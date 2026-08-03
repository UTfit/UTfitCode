import re
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("base_dir", help="Base directory, e.g. Summer26/NPWC_ICHEP")
args = parser.parse_args()

prefix = str(args.base_dir).replace("/", "_")
outfile = f"{prefix}-Lambda_summary.txt"
plotfile = f"{prefix}-Lambda_table.input"

interval_re = re.compile(r'\(\s*([-\deE.+]+)\s*,\s*([-\deE.+]+)\s*\)')

results = {}
max_lambda = {}

for i in range(1, 6):

    # Coefficients to look for in this file
    coeffs = {
        f"reC{i}_s",
        f"imC{i}_s",
        f"reC{i}_c",
        f"imC{i}_c",
        f"reC{i}_bd",
        f"imC{i}_bd",
        f"reC{i}_bs",
        f"imC{i}_bs",
    }
    
    filename = f"{args.base_dir}/C{i}/Observables/Statistics.txt"

    with open(filename) as f:
        lines = f.readlines()

    current_obs = None
    count = 0

    for j, line in enumerate(lines):

        m = re.search(r'Observable\s+"([^"]+)"', line)
        if m:
            current_obs = m.group(1)
            count = 0
            continue

        if current_obs not in coeffs:
            continue

        if "Smallest interval(s) containing at least" in line:
            count += 1

            if count == 2:
                m = interval_re.search(lines[j + 1])
                if m:
                    lower = float(m.group(1))
                    upper = float(m.group(2))

                    upperlimit = max(abs(lower), abs(upper))
                    lam = math.sqrt(1.0 / upperlimit) / 1000.

                    results[current_obs] = {
                        "upperlimit": upperlimit,
                        "lambda": lam,
                    }
                    # find the largest Lambda for this Ci
                    if i not in max_lambda or lam > max_lambda[i]["lambda"]:
                        max_lambda[i] = {
                            "observable": current_obs,
                            "upperlimit": upperlimit,
                            "lambda": lam,
                        }

print(f"{'Coefficient':12s} {'Upper limit':>15s} {'Lambda':>15s}")
print("-" * 45)
for coeff, lam in sorted(results.items()):
    upper = results[coeff]["upperlimit"]
    lam = results[coeff]["lambda"]
    print(f"{coeff:12s} {upper:15.6e} {lam:15.6e}")

for i in range(1, 6):
    entry = max_lambda[i]
    print(
        f"C{i:1d} "
        f"{entry['observable']:12s} "
        f"{entry['upperlimit']:15.6e} "
        f"{entry['lambda']:15.6e}"
    )


with open(outfile, "w") as out:

    out.write("All coefficients\n")
    out.write(f"{'Coefficient':12s} {'Upper limit':>15s} {'Lambda':>15s}\n")
    out.write("-" * 45 + "\n")

    for coeff in sorted(results):
        upper = results[coeff]["upperlimit"]
        lam = results[coeff]["lambda"]
        out.write(f"{coeff:12s} {upper:15.6e} {lam:15.6e}\n")

    out.write("\n")
    out.write("Maximum Lambda for each operator\n")
    out.write(f"{'C':>3s} {'Observable':12s} {'Upper limit':>15s} {'Lambda':>15s}\n")
    out.write("-" * 52 + "\n")

    for i in range(1, 6):
        entry = max_lambda[i]
        out.write(
            f"C{i:1d} "
            f"{entry['observable']:12s} "
            f"{entry['upperlimit']:15.6e} "
            f"{entry['lambda']:15.6e}\n"
        )

print(f"Results written to {outfile}")


labels = [
    ("Re(C_K)",  "re", "s"),
    ("Im(C_K)",  "im", "s"),
    ("Re(C_D)",  "re", "c"),
    ("Im(C_D)",  "im", "c"),
    ("Re(C_{B_d})", "re", "bd"),
    ("Im(C_{B_d})", "im", "bd"),
    ("Re(C_{B_s})", "re", "bs"),
    ("Im(C_{B_s})", "im", "bs")
]

with open(plotfile, "w") as out:

    out.write("scenario: Generic Flavour Structure\n")
    out.write("operators: $C_1$, $C_2$, $C_3$, $C_4$, $C_5$\n")

    for label, part, suffix in labels:

        values = []

        for i in range(1, 6):
            coeff = f"{part}C{i}_{suffix}"
            values.append(results[coeff]["lambda"])

        values_str = ", ".join(f"{v:.1f}" for v in values)

        out.write(f"$\\mathrm{{{label}}}$, {values_str}\n")

print(f"Table written to {plotfile}")
