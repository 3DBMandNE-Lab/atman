# Atman feature requests

No open requests. All filed items from the 2026-04-19 CSF cross-disease
batch shipped the same day.

## Notes for atman maintainers

- The CSF_CrossDisease manuscript is the first external consumer of atman
  v1.0.0 beyond the Dube NPX reproduction path. Any feature work here
  benefits other proteomics-analysis consumers of the canonical TSV
  contract.
- Contact: see `/Users/kevinjoseph/Cursor/CSF_CrossDisease/manuscript/`
  for the downstream analytic context each feature targets.
- Seed: all CSF manuscript analyses use `20260418` for direct
  compatibility with existing Python output; keeping that as the default
  seed in bootstrap/ICA/null commands simplifies cross-validation.
