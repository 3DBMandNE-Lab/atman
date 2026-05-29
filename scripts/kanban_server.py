#!/usr/bin/env python3
"""Kanban status server for biopsyID multi-agent orchestration.

See Docs/specs/2026-05-18-multi-agent-kanban-design.md for full design.
"""
from __future__ import annotations

import re
from typing import Optional

# ===========================================================================
# 1. Parser — read/write the kanban-status block in task markdown files.
# ===========================================================================

_BLOCK_RE = re.compile(
    r"<!--\s*kanban-status:start\s*-->(.*?)<!--\s*kanban-status:end\s*-->",
    re.DOTALL,
)

_FIELDS = ("status", "owner", "branch", "started", "updated")


def parse_status_block(md: str) -> Optional[dict]:
    """Return the parsed fields between the markers, or None if absent."""
    m = _BLOCK_RE.search(md)
    if not m:
        return None
    body = m.group(1)
    out = {f: None for f in _FIELDS}
    for line in body.splitlines():
        line = line.strip()
        for f in _FIELDS:
            key = f"**{f.capitalize()}:**"
            if line.startswith(key):
                out[f] = line[len(key):].strip()
                break
    return out


def write_status_block(md: str, fields: dict) -> str:
    """Rewrite the marker block contents. Raises ValueError if no block found."""
    m = _BLOCK_RE.search(md)
    if not m:
        raise ValueError("no kanban-status block found")
    lines = ["", *[
        f"**{f.capitalize()}:** {fields.get(f, '') or ''}"
        for f in _FIELDS
        if f in fields
    ], ""]
    new_block = (
        "<!-- kanban-status:start -->"
        + "\n".join(lines)
        + "<!-- kanban-status:end -->"
    )
    return md[:m.start()] + new_block + md[m.end():]


def make_claim_markdown(task_id: str, title: str = "",
                        owner: str = "orchestrator",
                        started: Optional[str] = None) -> str:
    """Generate canonical claim markdown for a new task.

    The orchestrator MUST use this rather than hand-typing the markdown, because
    fields written in non-canonical order create a merge-base whose layout
    differs from what `write_status_block` produces. When a subagent or the
    server rewrites the block later, git's recursive merger sees the whole
    block as changed on both sides and conservatively flags a conflict.

    Returns the full file content (caller writes it to tasks/<task_id>.md).
    """
    import datetime as _dt
    if started is None:
        started = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    heading = f"# {task_id}" + (f" — {title}" if title else "")
    template = (
        heading + "\n\n"
        "<!-- kanban-status:start -->\n"
        "<!-- kanban-status:end -->\n\n"
        "## Plan\n\n(filled in by subagent)\n\n"
        "## Notes\n\n(filled in by subagent)\n\n"
        "## Verification\n\n(filled in by subagent)\n\n"
        "## Hand-off\n\n(filled in by subagent)\n"
    )
    return write_status_block(template, {
        "status": "in_progress",
        "owner": owner,
        "branch": "(pending)",
        "started": started,
        "updated": started,
    })


# ===========================================================================
# 2. State builder — walk tasks/*.md on main, then read each task from the
#    right ref (branch if active, else main). See spec §5.3.
# ===========================================================================

import subprocess
from pathlib import Path

_TASK_FILENAME_RE = re.compile(r"^(TASK-\d+)\.md$")


def _git_quiet(args, cwd):
    """Run git, return (returncode, stdout, stderr); never raises."""
    res = subprocess.run(
        ["git", *args], cwd=cwd, capture_output=True, text=True,
    )
    return res.returncode, res.stdout, res.stderr


def _git_committer_iso(root, relpath):
    """Return ISO-8601 committer date of the last commit touching the file."""
    rc, out, _ = _git_quiet(["log", "-1", "--format=%cI", "--", relpath], cwd=root)
    if rc != 0:
        return None
    return out.strip() or None


def _branch_exists(root, branch):
    rc, _, _ = _git_quiet(["rev-parse", "--verify", f"refs/heads/{branch}"], cwd=root)
    return rc == 0


def _branch_is_merged_to_main(root, branch):
    rc, _, _ = _git_quiet(["merge-base", "--is-ancestor", branch, "main"], cwd=root)
    return rc == 0


def _read_task_md_from_branch(root, branch, task_id):
    """Return the markdown text from <branch>:tasks/TASK-XXX.md, or None."""
    rc, out, _ = _git_quiet(["show", f"{branch}:tasks/{task_id}.md"], cwd=root)
    if rc != 0:
        return None
    return out


_KANBAN_TASK_RE = re.compile(
    r'\{id:"(TASK-\d+)".*?depends_on:\[([^\]]*)\]',
    re.DOTALL,
)


def _parse_kanban_deps(root):
    """Parse kanban.html's task list and return {task_id: [dep_ids]}.

    Returns {} when kanban.html is absent or unreadable. Used by build_state
    to compute 'blocked' status for tasks that have no claim file.
    """
    p = Path(root) / "kanban.html"
    try:
        src = p.read_text()
    except OSError:
        return {}
    out = {}
    for m in _KANBAN_TASK_RE.finditer(src):
        tid = m.group(1)
        deps = re.findall(r'"(TASK-\d+)"', m.group(2))
        out[tid] = deps
    return out


def build_state(root):
    """Walk main's tasks/*.md, then read each from branch if appropriate.

    Additionally synthesise entries for any TASK-XXX listed in kanban.html
    but without a claim file: marked 'blocked' if any dependency is not yet
    done, else 'backlog'. This lets the dashboard's blocked column surface
    real dependency state instead of being empty.
    """
    root = Path(root)
    tasks_dir = root / "tasks"
    if not tasks_dir.is_dir():
        return {}
    state = {}
    for p in sorted(tasks_dir.iterdir()):
        m = _TASK_FILENAME_RE.match(p.name)
        if not m:
            continue
        task_id = m.group(1)

        # First pass: read main's view to discover the canonical Branch field.
        try:
            main_md = p.read_text()
        except OSError:
            continue
        main_fields = parse_status_block(main_md)
        canonical_branch = (main_fields or {}).get("branch")

        # Second pass: if branch is recorded, exists, and is not already
        # merged into main, read the status block from the branch's revision.
        fields = main_fields
        if (
            canonical_branch
            and canonical_branch != "(pending)"
            and _branch_exists(root, canonical_branch)
            and not _branch_is_merged_to_main(root, canonical_branch)
        ):
            branch_md = _read_task_md_from_branch(root, canonical_branch, task_id)
            if branch_md is not None:
                branch_fields = parse_status_block(branch_md)
                if branch_fields is not None:
                    fields = branch_fields

        if fields is None:
            state[task_id] = {
                "status": "unknown",
                "owner": None,
                "branch": canonical_branch,
                "claimed_at": None,
                "updated_at": _git_committer_iso(root, f"tasks/{p.name}"),
            }
            continue

        state[task_id] = {
            "status": fields["status"],
            "owner": fields["owner"],
            # `branch` returned to the UI is always main's canonical record,
            # never the branch's view of itself.
            "branch": canonical_branch,
            "claimed_at": fields["started"],
            "updated_at": fields["updated"] or _git_committer_iso(root, f"tasks/{p.name}"),
        }

    # Synthesise entries for tasks declared in kanban.html but without a
    # claim file. Status = 'blocked' if any dep is not yet done, else 'backlog'.
    deps_map = _parse_kanban_deps(root)
    done_ids = {tid for tid, entry in state.items() if entry.get("status") == "done"}
    for task_id, deps in deps_map.items():
        if task_id in state:
            continue
        missing = [d for d in deps if d not in done_ids]
        state[task_id] = {
            "status": "blocked" if missing else "backlog",
            "owner": None,
            "branch": None,
            "claimed_at": None,
            "updated_at": None,
            "missing_deps": missing or None,
        }
    return state


# ===========================================================================
# 3. Transitions — pure-function validator for status changes.
# ===========================================================================

# allowed_ui_transitions: (from_status, to_status) pairs the HTTP POST endpoint
# will accept. Other actors (orchestrator, subagent) have their own paths that
# do not go through the HTTP server.
_UI_ALLOWED = {
    ("review", "done"),
    ("review", "in_progress"),
    ("review", "blocked"),
    ("blocked", "in_progress"),
}

_VALID_STATUSES = {"backlog", "in_progress", "review", "blocked", "done"}


def validate_transition(current, target, actor):
    """Return None if allowed, else a failure code string."""
    if actor != "ui":
        return "invalid_transition"
    if target not in _VALID_STATUSES:
        return "invalid_transition"
    if (current, target) not in _UI_ALLOWED:
        return "invalid_transition"
    return None


def preflight_done(root, task_id):
    """Return None if the task is safe to promote to done, else (code, detail).

    Reads status FROM THE BRANCH if a branch is recorded — same read logic as
    build_state. This is critical: the subagent's `Status: review` lives on
    its task branch, not on main.
    """
    root = Path(root)
    task_path = root / "tasks" / f"{task_id}.md"
    if not task_path.exists():
        return ("not_found", f"{task_path} does not exist")

    # Get the canonical branch from main's record
    main_md = task_path.read_text()
    main_fields = parse_status_block(main_md) or {}
    branch = main_fields.get("branch")

    if not branch or branch == "(pending)":
        return ("branch_unrecorded", "**Branch:** still (pending); subagent never recorded it")
    rc, _, _ = _git_quiet(["rev-parse", "--verify", f"refs/heads/{branch}"], cwd=root)
    if rc != 0:
        return ("branch_missing", f"branch {branch!r} does not exist in this repo")

    # Now read the *current* status from the branch
    branch_md = _read_task_md_from_branch(root, branch, task_id)
    if branch_md is None:
        return ("not_found", f"branch {branch!r} has no tasks/{task_id}.md")
    branch_fields = parse_status_block(branch_md)
    if not branch_fields or branch_fields.get("status") != "review":
        cur = branch_fields.get("status") if branch_fields else None
        return ("invalid_transition",
                f"current status on branch is {cur!r}, expected 'review'")

    # main clean? Use --untracked-files=no so .gitignore'd files (events log)
    # don't poison the check; the .gitignore in Task 1 also keeps it out.
    rc, out, _ = _git_quiet(
        ["status", "--porcelain", "--untracked-files=no"], cwd=root,
    )
    if rc != 0 or out.strip():
        return ("dirty_main", "main working tree has uncommitted changes")
    return None


class BranchInUse(Exception):
    """Raised when git worktree add refuses because the branch is already
    checked out in another worktree (typically an unpruned Agent worktree)."""


def commit_to_branch(root, branch, repo_relpath, new_content, message):
    """Commit `new_content` to `repo_relpath` on `branch` via a transient
    worktree, without disturbing the main working tree.

    Raises BranchInUse if `git worktree add` fails because the branch is
    checked out elsewhere. Raises RuntimeError on other git failures.
    """
    import tempfile
    root = Path(root)
    # Clear stale worktree refs (where the dir was deleted manually).
    _git_quiet(["worktree", "prune"], cwd=root)

    tmpdir = tempfile.mkdtemp(prefix="kanban-srv-")
    try:
        rc, _, err = _git_quiet(["worktree", "add", tmpdir, branch], cwd=root)
        if rc != 0:
            # Git says "is already checked out at ..." if the branch is taken.
            if "already checked out" in err.lower() or "already used by worktree" in err.lower():
                raise BranchInUse(err.strip())
            raise RuntimeError(f"git worktree add failed: {err.strip()}")

        target = Path(tmpdir) / repo_relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(new_content)

        rc, _, err = _git_quiet(["add", repo_relpath], cwd=tmpdir)
        if rc != 0:
            raise RuntimeError(f"git add failed: {err.strip()}")
        rc, _, err = _git_quiet(["commit", "-m", message], cwd=tmpdir)
        if rc != 0:
            raise RuntimeError(f"git commit failed: {err.strip()}")
    finally:
        # Best-effort cleanup; ignore errors so a failed worktree-add does not
        # mask the original exception.
        _git_quiet(["worktree", "remove", "--force", tmpdir], cwd=root)


# ===========================================================================
# 4. Events log — append-only TSV operator log (not regulatory; see spec §7).
# ===========================================================================

import datetime as _dt


def _escape_tsv(s):
    """Replace tab and newline so a note cannot break TSV row/column layout."""
    if s is None:
        return ""
    return str(s).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def append_event(root, task_id, from_status, to_status, actor, note):
    log = Path(root) / "tasks" / ".events.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    ts = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    row = "\t".join(_escape_tsv(x) for x in [ts, task_id, from_status, to_status, actor, note])
    with log.open("a") as f:
        f.write(row + "\n")


# ===========================================================================
# 5. HTTP — request handler factory and main entry.
# ===========================================================================

import json
from http.server import BaseHTTPRequestHandler


def make_handler(root):
    """Return a request-handler class bound to a specific repo root.

    Returning a class (not a closure-handler) keeps the BaseHTTPRequestHandler
    contract straightforward — http.server instantiates it per request.
    """
    root_path = Path(root).resolve()

    class KanbanHandler(BaseHTTPRequestHandler):
        # quieter logs
        def log_message(self, fmt, *args):
            return

        def _send_json(self, status_code, payload):
            body = json.dumps(payload).encode()
            self.send_response(status_code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _send_html(self, status_code, body_str):
            body = body_str.encode()
            self.send_response(status_code)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path == "/" or self.path == "/index.html":
                kanban = root_path / "kanban.html"
                if not kanban.exists():
                    self._send_json(500, {"error": "kanban_html_missing"})
                    return
                self._send_html(200, kanban.read_text())
                return
            if self.path == "/state.json":
                self._send_json(200, build_state(root_path))
                return
            self._send_json(404, {"error": "not_found", "detail": self.path})

        def _read_body(self):
            try:
                length = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                length = 0
            raw = self.rfile.read(length) if length else b""
            try:
                return json.loads(raw or b"{}")
            except json.JSONDecodeError:
                return None

        def _current_branch_status(self, task_id, main_fields):
            """Return (current_status, branch_or_None) using the same read
            logic as build_state — branch if active, else main."""
            branch = (main_fields or {}).get("branch")
            if (
                branch and branch != "(pending)"
                and _branch_exists(root_path, branch)
                and not _branch_is_merged_to_main(root_path, branch)
            ):
                branch_md = _read_task_md_from_branch(root_path, branch, task_id)
                if branch_md is not None:
                    bf = parse_status_block(branch_md)
                    if bf is not None:
                        return bf.get("status"), branch
            return (main_fields or {}).get("status"), branch

        def do_POST(self):
            m = re.match(r"^/tasks/(TASK-\d+)/status/?$", self.path)
            if not m:
                self._send_json(404, {"error": "not_found", "detail": self.path})
                return
            task_id = m.group(1)

            body = self._read_body()
            if body is None or not isinstance(body, dict) or "status" not in body:
                self._send_json(400, {"error": "bad_request",
                                       "detail": "body must be JSON with a 'status' field",
                                       "task_id": task_id})
                return

            target = body["status"]
            note = body.get("note", "") or ""

            task_path = root_path / "tasks" / f"{task_id}.md"
            if not task_path.exists():
                self._send_json(404, {"error": "not_found", "task_id": task_id})
                return

            main_md = task_path.read_text()
            main_fields = parse_status_block(main_md) or {}
            current, branch = self._current_branch_status(task_id, main_fields)

            if target == "done":
                pf = preflight_done(root_path, task_id)
                if pf is not None:
                    pcode, detail = pf
                    http_code = {
                        "not_found": 404,
                        "invalid_transition": 400,
                        "branch_unrecorded": 409,
                        "branch_missing": 409,
                        "dirty_main": 409,
                    }.get(pcode, 400)
                    self._send_json(http_code, {"error": pcode, "detail": detail,
                                                "task_id": task_id})
                    return

                # Already-merged short-circuit: if the branch is an ancestor of
                # main, the merge step is a no-op; just rewrite status to done.
                rc, _, _ = _git_quiet(
                    ["merge-base", "--is-ancestor", branch, "main"], cwd=root_path,
                )
                if rc != 0:
                    # Real merge needed.
                    rc, _, err = _git_quiet(
                        ["merge", "--no-ff", branch,
                         "-m", f"kanban: {task_id} → done (merge {branch})"],
                        cwd=root_path,
                    )
                    if rc != 0:
                        _git_quiet(["merge", "--abort"], cwd=root_path)
                        self._send_json(422, {"error": "merge_conflict",
                                               "detail": err.strip()[:1000],
                                               "task_id": task_id})
                        return

                # Capture the (possibly short-circuited) HEAD sha.
                rc, sha, _ = _git_quiet(["rev-parse", "HEAD"], cwd=root_path)
                merged_commit = sha.strip() if rc == 0 else None

                # Rewrite the status block on main to done. The merge brought in
                # the branch's markdown (Status: review); we now overwrite to
                # Status: done and commit on main.
                now = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
                md_after = task_path.read_text()
                fields_after = parse_status_block(md_after) or main_fields
                new_fields = dict(fields_after)
                new_fields["status"] = "done"
                new_fields["updated"] = now
                new_md = write_status_block(md_after, new_fields)
                task_path.write_text(new_md)
                _git_quiet(["add", f"tasks/{task_id}.md"], cwd=root_path)
                rc, _, err = _git_quiet(
                    ["commit", "-m", f"kanban: {task_id} → done"], cwd=root_path,
                )
                if rc != 0:
                    if "nothing to commit" not in err.lower():
                        self._send_json(500, {"error": "commit_failed",
                                               "detail": err.strip(),
                                               "task_id": task_id})
                        return

                # Best-effort worktree prune.
                _git_quiet(["worktree", "prune"], cwd=root_path)

                append_event(root_path, task_id, "review", "done", "ui", note)
                self._send_json(200, {"task_id": task_id, "status": "done",
                                       "branch": branch,
                                       "merged_commit": merged_commit})
                return

            code = validate_transition(current, target, "ui")
            if code is not None:
                self._send_json(400, {"error": code, "task_id": task_id,
                                       "detail": f"{current!r} → {target!r} not permitted"})
                return

            # Non-merge transition lands on the branch (where the read path looks).
            # If the branch is gone (e.g. blocked → in_progress after manual cleanup),
            # fall through to writing on main.
            now = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
            write_to_main = not (branch and branch != "(pending)" and _branch_exists(root_path, branch))

            if write_to_main:
                new_fields = dict(main_fields)
                new_fields["status"] = target
                new_fields["updated"] = now
                new_md = write_status_block(main_md, new_fields)
                task_path.write_text(new_md)
                _git_quiet(["add", f"tasks/{task_id}.md"], cwd=root_path)
                rc, _, err = _git_quiet(
                    ["commit", "-m", f"kanban: {task_id} → {target} (on main, branch gone)"],
                    cwd=root_path,
                )
                if rc != 0:
                    self._send_json(500, {"error": "commit_failed", "task_id": task_id,
                                           "detail": err.strip()})
                    return
            else:
                # Read the current branch state, edit, commit via transient worktree.
                branch_md = _read_task_md_from_branch(root_path, branch, task_id) or main_md
                branch_fields = parse_status_block(branch_md) or main_fields
                new_fields = dict(branch_fields)
                new_fields["status"] = target
                new_fields["updated"] = now
                new_md = write_status_block(branch_md, new_fields)
                try:
                    commit_to_branch(root_path, branch, f"tasks/{task_id}.md", new_md,
                                     f"kanban: {task_id} → {target}")
                except BranchInUse as e:
                    self._send_json(409, {"error": "branch_in_use",
                                           "task_id": task_id,
                                           "detail": str(e)[:500]})
                    return
                except RuntimeError as e:
                    self._send_json(500, {"error": "commit_failed",
                                           "task_id": task_id,
                                           "detail": str(e)[:500]})
                    return

            append_event(root_path, task_id, current, target, "ui", note)
            self._send_json(200, {"task_id": task_id, "status": target})

    return KanbanHandler


def main(argv=None):
    import argparse, os
    p = argparse.ArgumentParser(description="biopsyID kanban status server")
    p.add_argument("--port", type=int, default=7878)
    p.add_argument("--root", default=os.getcwd())
    args = p.parse_args(argv)

    from http.server import HTTPServer
    handler = make_handler(args.root)
    server = HTTPServer(("127.0.0.1", args.port), handler)
    print(f"kanban-server: listening on http://127.0.0.1:{args.port} (root={args.root})",
          flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("kanban-server: shutting down")
        server.shutdown()


if __name__ == "__main__":
    main()
