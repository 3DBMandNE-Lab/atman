#!/usr/bin/env python3
"""Generate inputs and a reference table for tests/axes_contrast.rs.

Synthetic two-cohort score table; reference statistics from numpy/scipy/
statsmodels using the manuscript's definitions (pooled-SD Cohen d, Hedges
large-sample d CI, Welch t, Mann-Whitney AUC, OLS with z(age) + sex).
"""
from pathlib import Path
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats

OUT = Path(__file__).parent / "axes_contrast"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(7)
rows = []
for cohort, n, case_label, ctrl_label in [("c1", 14, "Case", "Ctrl"), ("c2", 12, "Dis", "Ref")]:
    for i in range(n):
        is_case = i < n // 2
        rows.append(dict(sample_id=f"{cohort}_{i:02d}", cohort=cohort, condition=case_label if is_case else ctrl_label,
                         is_control=0 if is_case else 1, age=float(rng.integers(30, 80)), sex=rng.choice(["F", "M"]),
                         site=["Berlin", "Kiel", "Sweden"][i % 3],
                         s1=rng.normal(0.8 if is_case else 0.0, 1.0), s2=rng.normal(0.0, 1.0)))
df = pd.DataFrame(rows)
df.loc[3, "age"] = np.nan          # one missing covariate in c1
df.loc[20, "s2"] = np.nan          # one missing score in c2
df[["sample_id", "cohort", "condition", "is_control", "s1", "s2"]].to_csv(OUT / "scores.tsv", sep="\t", index=False, na_rep="")
for cohort in ["c1", "c2"]:
    sub = df[df.cohort == cohort]
    d = OUT / cohort; d.mkdir(exist_ok=True)
    pd.DataFrame(dict(sample_id=sub.sample_id, subject_id=sub.sample_id, condition=sub.condition, is_control=sub.is_control,
                      sample_type="bio", ingest_order=range(1, len(sub) + 1), age=sub.age, sex=sub.sex, site=sub.site)).to_csv(d / "samples.tsv", sep="\t", index=False, na_rep="")
pd.DataFrame([dict(label="c1_case", cohort="c1", case="Case", control="Ctrl", family="fam"),
              dict(label="c2_case", cohort="c2", case="Dis", control="*", family="fam")]).to_csv(OUT / "contrasts.tsv", sep="\t", index=False)


def d_ci(a, b):
    n1, n2 = len(a), len(b)
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    d = (a.mean() - b.mean()) / sp
    se = np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    return d, d - 1.96 * se, d + 1.96 * se


ref = []
for label, cohort, case in [("c1_case", "c1", "Case"), ("c2_case", "c2", "Dis")]:
    sub = df[df.cohort == cohort].copy(); sub["case"] = (sub.condition == case).astype(float)
    sub["sex_m"] = (sub.sex == "M").astype(float)
    for score in ["s1", "s2"]:
        s = sub.dropna(subset=[score])
        a, b = s[s.case == 1][score].values, s[s.case == 0][score].values
        d, lo, hi = d_ci(a, b); t, p = stats.ttest_ind(a, b, equal_var=False)
        auc = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
        for design, terms in [("~ case", ["case"]), ("~ case + z(age) + sex", ["case", "age_s", "sex_m"])]:
            cc = s.dropna(subset=["age"]) if "age_s" in terms else s
            cc = cc.copy(); cc["age_s"] = (cc.age - cc.age.mean()) / cc.age.std(ddof=1)
            fit = sm.OLS(cc[score].astype(float), sm.add_constant(cc[terms].astype(float))).fit()
            ci = fit.conf_int().loc["case"]
            ref.append(dict(label=label, score=score, design=design, n_case=len(a), n_control=len(b), mean_diff=a.mean() - b.mean(),
                            cohen_d=d, d_ci_lo=lo, d_ci_hi=hi, welch_t=t, welch_p=p, auc=auc, n_fit=int(fit.nobs),
                            beta=fit.params["case"], se=fit.bse["case"], ci_lo=ci[0], ci_hi=ci[1], p=fit.pvalues["case"],
                            beta_over_resid_sd=fit.params["case"] / np.sqrt(fit.scale)))
pd.DataFrame(ref).to_csv(OUT / "reference.tsv", sep="\t", index=False, float_format="%.10f")

# Site sensitivity: c1_case, s1, ~ case + z(age) + sex + site ; joint F-test of the site dummies.
sub = df[df.cohort == "c1"].dropna(subset=["age"]).copy(); sub["case"] = (sub.condition == "Case").astype(float)
sub["sex_m"] = (sub.sex == "M").astype(float); sub["age_s"] = (sub.age - sub.age.mean()) / sub.age.std(ddof=1)
X = pd.DataFrame(dict(const=1.0, case=sub.case, age_s=sub.age_s, sex_m=sub.sex_m,
                      siteKiel=(sub.site == "Kiel").astype(float), siteSweden=(sub.site == "Sweden").astype(float)))
fit = sm.OLS(sub.s1.astype(float), X).fit()
ft = fit.f_test(np.array([[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]]))
pd.DataFrame([dict(beta_case=fit.params["case"], p_case=fit.pvalues["case"], f=float(np.squeeze(ft.fvalue)), df_num=float(ft.df_num), df_den=float(ft.df_denom), p_f=float(np.squeeze(ft.pvalue)), n_fit=int(fit.nobs))]).to_csv(OUT / "reference_site.tsv", sep="\t", index=False, float_format="%.10f")
print("wrote", OUT)
