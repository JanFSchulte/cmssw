"""Simple, transparent variant x eta-region comparison table.

No automated "winner" selection or optimization -- this only tabulates
response bias, resolution, and fake rate (at a fixed pt/PU working point) per
variant so a physicist can make the final call. Keep it that way; this is a
calibration aid, not a tuning algorithm.
"""

import numpy as np

from . import response as response_mod


def summarize_variant(table_by_variant, jets_by_variant, nTrueInt_by_variant,
                       genjet_pt_min=30.0, pt_cut=30.0, eta_max=2.5):
    """table_by_variant: {variant_name: response_table(...)} (see response.py)
    jets_by_variant / nTrueInt_by_variant: raw per-variant jets + PU arrays,
    needed separately for fake_rate_vs_pu (which works on unmatched jets too).

    Returns a list of dict rows: variant, response_bias (median response - 1,
    averaged over genjet_pt>genjet_pt_min), resolution (averaged IQR/2),
    fake_rate (at pt_cut).
    """
    rows = []
    for variant, table in table_by_variant.items():
        sel = table["genjet_pt"] > genjet_pt_min
        if not np.any(sel):
            continue
        bias = np.median(table["response"][sel]) - 1.0
        # resolution: IQR/2 of response in the same selection
        q16, q84 = np.percentile(table["response"][sel], [16, 84])
        resolution = 0.5 * (q84 - q16)

        fake_rate = float("nan")
        if variant in jets_by_variant and nTrueInt_by_variant.get(variant) is not None:
            _, rates, counts = response_mod.fake_rate_vs_pu(
                jets_by_variant[variant], nTrueInt_by_variant[variant],
                pt_cut=pt_cut, eta_max=eta_max,
            )
            if len(rates):
                fake_rate = float(np.average(rates, weights=counts))

        rows.append({
            "variant": variant,
            "response_bias": bias,
            "resolution": resolution,
            "fake_rate": fake_rate,
        })
    return rows


def format_markdown_table(rows):
    header = "| variant | response bias | resolution (IQR/2) | fake rate @ fixed pt |"
    sep = "|---|---|---|---|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            "| %s | %+.3f | %.3f | %.3f |"
            % (r["variant"], r["response_bias"], r["resolution"], r["fake_rate"])
        )
    return "\n".join(lines)


def write_csv(rows, path):
    import csv
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["variant", "response_bias", "resolution", "fake_rate"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
