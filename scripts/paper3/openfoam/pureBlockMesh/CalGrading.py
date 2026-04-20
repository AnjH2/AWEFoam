import math

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
        nr_=self.calculate('calculate_nr_from_cratio_widthA_length', cratio_, widthA_,length_)[0]
        ratio_ = self.calculate('calculate_ratio_from_nr_cratio', nr_, cratio_)
        return nr_,ratio_
# Example usage
if __name__ == "__main__":
    widthDi=0.063;
    widthCh=0.126;
    cratio=1.2;
    tCh=3;
    tEl=1.8;
    tDi=0.85;
    all_funcs = AllFunctions()
    known_values = {
        "cratio": 1.2,
    }

   
    
    print("known_values:", known_values)
    nChannel,rChannel = all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', known_values['cratio'], widthCh,tCh/2)
    print("Channel")
    print(nChannel*2)
    print(rChannel)
    nDia,rDia = all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', known_values['cratio'], widthDi,tDi/2)
    print("Dia")
    print(nDia*2)
    print(rDia)
    nEl,rEl = all_funcs.calculate('calculate_nr_and_ratio_from_cratio_widthA_length', known_values['cratio'], widthCh*2,tEl/2)
    print("el")
    print(nEl*2)
    print(rEl)


