#!/usr/bin/env python3
import math
import os


class AllFunctions:
    def __init__(self, precision=10, relTol=1e-5, rMax=1e5):
        self.precision = precision
        self.relTol = relTol
        self.rMax = rMax

    def sum_cells(self, x0, r, n):
        if n is None or n < 1:
            return 0
        else:
            length = 0
            x = x0
            for i in range(1, n + 1):
                length += x
                x *= r
            return length

    def root_by_bisection(self, f, x1, x2):
        max_steps = 500
        f1 = f(x1)
        f2 = f(x2)
        if f1 * f2 >= 0:
            print(f"{x1} and {x2} do not bracket the root ({f1}, {f2})")
            return None, "Internal error: Wrong start values for root finding"

        steps = 0
        while steps < max_steps:
            steps += 1
            x_mid = 0.5 * (x1 + x2)
            if (x2 - x1) < self.relTol:
                return x_mid
            f_mid = f(x_mid)
            if f1 * f_mid < 0:
                x2 = x_mid
                f2 = f_mid
            else:
                x1 = x_mid
                f1 = f_mid

        return None, "Internal error: root finding did not converge"

    def round_val(self, val):
        return [round(val), float(format(val, f".{self.precision}g"))]

    def calculate_ratio_from_nr_cratio(self, nr, cratio):
        if nr > 1:
            return math.pow(cratio, nr - 1)
        else:
            return [cratio, "Number of cells must be >1"]

    def calculate_cratio_from_nr_ratio(self, nr, ratio):
        if nr > 1:
            return math.pow(ratio, 1 / (nr - 1))
        else:
            return [ratio, "Number of cells must be >1"]

    def calculate_ratio_from_widthA_widthB(self, widthA, widthB):
        return widthB / widthA

    def calculate_widthB_from_widthA_ratio(self, widthA, ratio):
        return widthA * ratio

    def calculate_widthA_from_widthB_ratio(self, widthB, ratio):
        return widthB / ratio

    def calculate_length_from_widthA_cratio_nr(self, widthA, cratio, nr):
        return self.sum_cells(widthA, cratio, nr)

    def calculate_length_from_widthB_cratio_nr(self, widthB, cratio, nr):
        return self.sum_cells(widthB, 1 / cratio, nr)

    def calculate_nr_from_ratio_cratio(self, ratio, cratio):
        if abs(cratio - 1) > 1e-5:
            return self.round_val(math.log(ratio) / math.log(cratio) + 1)
        else:
            return [None, "Cell to cell ratio must not be 1"]

    def calculate_widthA_from_nr_cratio_length(self, nr, cratio, length):
        if abs(cratio - 1) > 1e-5:
            return length * (1 - cratio) / (1 - math.pow(cratio, nr))
        else:
            return length / nr

    def calculate_nr_from_cratio_widthA_length(self, cratio, widthA, length):
        if abs(cratio - 1) > self.relTol:
            return self.round_val(math.log(1 - length / widthA * (1 - cratio)) / math.log(cratio))
        else:
            return self.round_val(length / widthA)

    def calculate_nr_from_cratio_widthB_length(self, cratio, widthB, length):
        if abs(cratio - 1) > self.relTol:
            return self.round_val(math.log(1 / (1 + length / widthB * (1 - cratio) / cratio)) / math.log(cratio))
        else:
            return self.round_val(length / widthB)

    def calculate_cratio_from_nr_widthB_length(self, nr, widthB, length):
        if abs(nr * widthB - length) / length < self.relTol:
            return 1
        else:
            if nr * widthB > length:
                cMax = math.pow(self.rMax, 1 / (nr - 1))
                cMin = math.pow(1 + self.relTol, 1 / (nr - 1))
            else:
                cMax = math.pow(1 - self.relTol, 1 / (nr - 1))
                cMin = math.pow(1 / self.rMax, 1 / (nr - 1))

            return self.root_by_bisection(
                lambda c: (1 / math.pow(c, nr - 1)) * (1 - math.pow(c, nr)) / (1 - c) - length / widthB,
                cMin,
                cMax
            )

    def calculate_cratio_from_nr_widthA_length(self, nr, widthA, length):
        if abs(nr * widthA - length) / length < self.relTol:
            return 1
        else:
            if nr * widthA < length:
                cMax = math.pow(self.rMax, 1 / (nr - 1))
                cMin = math.pow(1 + self.relTol, 1 / (nr - 1))
            else:
                cMax = math.pow(1 - self.relTol, 1 / (nr - 1))
                cMin = math.pow(1 / self.rMax, 1 / (nr - 1))

            return self.root_by_bisection(
                lambda c: (1 - math.pow(c, nr)) / (1 - c) - length / widthA,
                cMin, cMax
            )

    def calculate_nr_from_ratio_widthA_length(self, ratio, widthA, length):
        if abs(ratio - 1) < self.relTol:
            return self.round_val(length / widthA)
        else:
            return self.round_val(
                self.root_by_bisection(
                    lambda n: (1 - math.pow(ratio, n / (n - 1))) / (1 - math.pow(ratio, 1 / (n - 1))) - length / widthA,
                    0, length / widthA
                )
            )

    def calculate(self, function_name, *args):
        if hasattr(self, function_name):
            return getattr(self, function_name)(*args)
        else:
            raise ValueError(f"Function {function_name} not found")

    def calculate_nr_and_ratio_from_cratio_widthA_length(self, cratio_, widthA_, length_):
        nr_ = self.calculate('calculate_nr_from_cratio_widthA_length', cratio_, widthA_, length_)[0]
        ratio_ = self.calculate('calculate_ratio_from_nr_cratio', nr_, cratio_)
        return nr_, ratio_


def grading_string_symmetric(ratio_half):
    return (
        "(\n"
        f"    (0.5 0.5 {ratio_half:.12g})\n"
        f"    (0.5 0.5 {1.0/ratio_half:.12g})\n"
        ")"
    )


def main():
    all_funcs = AllFunctions()

    # -------------------------------------------------
    # USER INPUTS
    # lengths in mm
    # -------------------------------------------------
    tCh = 3.0
    tEl = 1.8
    tDi = 0.85

    # target first-cell widths for the half-blocks, in mm
    # these match your example structure
    widthDi = 0.063
    widthCh = 0.063
    widthEl = widthCh

    # cell-to-cell ratio in the half-blocks
    cratio = 1.2

    # output file
    out_file = os.path.join(".","system", "meshZgrading.auto")

    # -------------------------------------------------
    # CALCULATIONS
    # -------------------------------------------------

    # z01 / z45 : electrode-style symmetric block
    nEl_half, rEl = all_funcs.calculate(
        'calculate_nr_and_ratio_from_cratio_widthA_length',
        cratio, widthEl, tEl / 2.0
    )
    n_z12 = int(nEl_half * 2)
    n_z34 = int(nEl_half * 2)

    # z12 / z34 : channel-style symmetric block
    nCh_half, rCh = all_funcs.calculate(
        'calculate_nr_and_ratio_from_cratio_widthA_length',
        cratio, widthCh, tCh / 2.0
    )
    n_z01 = int(nCh_half * 2)
    n_z45 = int(nCh_half * 2)

    # z23 : diaphragm-style symmetric block
    nDi_half, rDi = all_funcs.calculate(
        'calculate_nr_and_ratio_from_cratio_widthA_length',
        cratio, widthDi, tDi / 2.0
    )
    n_z23 = int(nDi_half * 2)

    gradZ01 = grading_string_symmetric(rEl)
    gradZ12 = grading_string_symmetric(rCh)
    gradZ23 = grading_string_symmetric(rDi)
    gradZ34 = grading_string_symmetric(rCh)
    gradZ45 = grading_string_symmetric(rEl)

    # reconstructed interface-to-center lengths for checking
    Lhalf01 = all_funcs.calculate('calculate_length_from_widthA_cratio_nr', widthEl, cratio, nEl_half)
    Lhalf12 = all_funcs.calculate('calculate_length_from_widthA_cratio_nr', widthCh, cratio, nCh_half)
    Lhalf23 = all_funcs.calculate('calculate_length_from_widthA_cratio_nr', widthDi, cratio, nDi_half)

    # -------------------------------------------------
    # WRITE FILE
    # -------------------------------------------------
    lines = []
    lines.append("// auto-generated by make_zgrading.py")
    lines.append("// all lengths and widths below are in mm")
    lines.append("")
    
    lines.append(f"tEl {tEl:.12g};")
    lines.append(f"tCh {tCh:.12g};")
    lines.append(f"tDi {tDi:.12g};")

    lines.append("// input lengths")
    lines.append(f"L_z01 {tEl:.12g};")
    lines.append(f"L_z12 {tCh:.12g};")
    lines.append(f"L_z23 {tDi:.12g};")
    lines.append(f"L_z34 {tCh:.12g};")
    lines.append(f"L_z45 {tEl:.12g};")
    lines.append("")

    lines.append("// input first-cell widths for half-block calculations")
    lines.append(f"widthEl {widthEl:.12g};")
    lines.append(f"widthCh {widthCh:.12g};")
    lines.append(f"widthDi {widthDi:.12g};")
    lines.append(f"cratioZ {cratio:.12g};")
    lines.append("")

    lines.append("// solved cell counts")
    lines.append(f"n_z01 {n_z01};")
    lines.append(f"n_z12 {n_z12};")
    lines.append(f"n_z23 {n_z23};")
    lines.append(f"n_z34 {n_z34};")
    lines.append(f"n_z45 {n_z45};")
    lines.append("")

    lines.append("// solved grading ratios for each half-block")
    lines.append(f"rEl {rEl:.12g};")
    lines.append(f"rCh {rCh:.12g};")
    lines.append(f"rDi {rDi:.12g};")
    lines.append("")

    lines.append("// reconstructed half lengths for checking")
    lines.append(f"Lhalf01_check {Lhalf01:.12g};")
    lines.append(f"Lhalf12_check {Lhalf12:.12g};")
    lines.append(f"Lhalf23_check {Lhalf23:.12g};")
    lines.append("")

    lines.append("gradZ01")
    lines.append(gradZ01 + ";")
    lines.append("")

    lines.append("gradZ12")
    lines.append(gradZ12 + ";")
    lines.append("")

    lines.append("gradZ23")
    lines.append(gradZ23 + ";")
    lines.append("")

    lines.append("gradZ34")
    lines.append(gradZ34 + ";")
    lines.append("")

    lines.append("gradZ45")
    lines.append(gradZ45 + ";")
    lines.append("")

    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"Wrote {out_file}")
    print()
    print("Summary")
    print("-------")
    print(f"z01/z45: nHalf={nEl_half}, n={n_z01}, ratioHalf={rEl}")
    print(f"z12/z34: nHalf={nCh_half}, n={n_z12}, ratioHalf={rCh}")
    print(f"z23    : nHalf={nDi_half}, n={n_z23}, ratioHalf={rDi}")
    print()
    print("Half-length checks")
    print("------------------")
    print(f"z01/z45 target={tEl/2.0}  check={Lhalf01}")
    print(f"z12/z34 target={tCh/2.0}  check={Lhalf12}")
    print(f"z23     target={tDi/2.0}  check={Lhalf23}")


if __name__ == "__main__":
    main()
