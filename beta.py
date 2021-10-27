import re
import copy
from lincomb import *

class ANSSData(object):
    def __init__(self, file_prefix="data/185"):
        self.beta_mapping = {
                '_BPBocSS_theta2.txt': 'beta1',
                '_BPBocSS_theta3.txt': 'beta2',
                '_BPBocSS_theta4.txt': 'beta33',
                '_BPBocSS_theta5.txt': 'beta4',
                '_BPBocSS_theta6.txt': 'beta5',
                '_BPBocSS_theta7.txt': 'beta63',
                }

        # E.g. the (relative) pathname file_prefix + "_BPBocSS_table.txt" should exist
        self.file_prefix = file_prefix
        # bockstein is a dict where key: value means d_r(key) = value (both strings, where "" means 0)
        # bottom_cells is a list of permanent cycles (as strings) in the Bockstein sseq
        self.bockstein, self.bottom_cells = self.make_bockstein()
        # boc_a0table is a dict where key: value means key * 3 = value as
        # permanent cycles in the Bockstein spectral sequence (so integral
        # classes, but Bockstein names)
        self.boc_a0table = self.parse_mult(self.file_prefix + "_BPBocSS_a0.txt")
        self.boc_h0table = self.parse_mult(self.file_prefix + "_BPBocSS_h0.txt")
        # The rest of these are straightforward multiplication tables
        self.a0table = self.parse_mult(self.file_prefix + "_BPAANSS_a0.txt")
        self.h0table = self.parse_mult(self.file_prefix + "_BPAANSS_h0.txt")
        self.beta_table = { v : self.parse_mult(self.file_prefix + k)
                for k,v in self.beta_mapping.items() }
        # deg_table is a dict where keys are names in E_2(S) and values are
        # dicts with keys 's' (stem), 'f' (ANSS filtration), 'nov' (Novikov
        # degree, i.e. Adams filtration = f + nov)
        self.deg_table = self.make_deg_table()
        # B2A is a dict where key:value means value is the algNSS name
        # corresponding to the Bockstein element key
        # (key is a string, value is a LinComb object)
        self.B2A = self.make_B2A()
        # B2A mapping reversed (keys are strings, values are LinComb)
        self.B2A_inv = self.make_B2A_inv()
        self.boc_a0div = self.make_boc_a0div()


    def make_bockstein(self):
        # BocSS_table contains no linear combinations, so bockstein is a
        # strings:strings mapping. This makes it easy to check whether an
        # element is in the image of the bockstein mapping.
        fh = open(self.file_prefix + "_BPBocSS_table.txt", "r")
        bockstein = {}
        bottom_cells = []
        for line in fh:
          if re.search("<-", line):
              target, _, src, dlen, _ = re.split("	", line)
              if not re.search("d0", dlen):
                  bockstein[src] = target
          else:
              elt, _ = re.split("	", line)
              bottom_cells.append(elt)
              bockstein[elt] = ""
        fh.close()
        return (bockstein, bottom_cells)

    def make_boc_a0table(self):
        # different from parse_mult because terms involved in d0's need to be removed
        fh = open(self.file_prefix + "_BPBocSS_a0.txt", "r")
        a0table = {}
        for line in fh:
            src, target = re.split("	->	", line)
            target = re.sub("\n", "", target)
            target_lc = LinComb.zero()
            for t in re.split("\+", target):
                if t in self.bottom_cells or t in self.bockstein.values():
                    target_lc.add_inplace(LinComb({t:1}))
            a0table[src] = target_lc
        fh.close()
        return a0table

    # Translation between different names for elements in E_2(S/p)
    def make_B2A(self):
        B2A = {}
        for line in open(self.file_prefix + "_BPB2A_table.txt"):
            line = re.sub("\n", "", line)
            pair = re.split("	->	", line)
            newval = LinComb.zero()
            terms = re.split("\+", pair[1])
            for t in terms:
                if t != "o":
                    newval.add_inplace(LinComb({t:1}))
            B2A[pair[0]] = newval
        return B2A

    def make_deg_table(self):
        filename = self.file_prefix + "_BPAANSS_table.txt"
        fh = open(filename, "r")
        rtn = {}
        for line in fh:
            if re.search("<-", line):
                continue
            items = re.split("\t\|", line)
            deg_pair = re.sub("[^0-9,]", "", items[-1])
            deg_pair = re.split(",", deg_pair)
            stem = int(deg_pair[0])
            adams_filt = int(deg_pair[1])
            # An element named v[a-b] will have ANSS filtration = a
            anss_filt = int(re.search("^[0-9]*", re.sub("[v^0-9]*\[", "", items[0])).group())
            nov = adams_filt - anss_filt
            name = items[0]
            rtn[name] = {'s': stem, 'f': anss_filt, 'nov': nov}
        return rtn

    def make_B2A_inv(self):
        max_f = max([r['f'] + r['nov'] for i,r in self.deg_table.items()])
        a0div = [gen for v in self.a0table.values() for gen in v.elts.keys()] # a0-divisible generators (excludes linear combinations)
        # Exclude a0-divisible elements since those go to zero along S --> S/3.
        # Also exclude values near the boundary of computation (since those seem to have issues).
        B2A_inv = self.invert(self.B2A,
                lambda x : (self.deg_table[x]['f'] + self.deg_table[x]['nov'] >= max_f
                    or x in a0div))
        # We removed the a0-divisible summands to make B2A invertible.
        # Actually, B2A_inv applied to all these elements should be zero.
        for i in self.deg_table.keys():
            keys = list(B2A_inv.keys())
            if not i in keys:
                B2A_inv[i] = LinComb.zero()
        return B2A_inv

    def make_boc_a0div(self):
        # Exclude terms involved in a bockstein d0 since those are gone by E_1
        return self.invert(self.boc_a0table,
                lambda x : not(x in self.bockstein.keys() or x in self.bockstein.values()))

    def invert(self, orig, exclude_fn=None):
        # Ad hoc approach to inverting the dict orig (containing string: LinComb pairs)
        # A term t in the values will be removed if exclude_fn(t) == True
        inv = {}
        for k,v in orig.items():
            targets = list(v.elts.keys()) # a list of generators involved in the linear combination orig[k]
            if len(targets) == 0:
                continue # This case shouldn't happen, but sometimes it does because k is close to the edge of the computation.
            elif len(targets) == 1:
                if targets[0] in inv.keys():
                    raise Exception("Failed to invert dict; %s is hit multiple times" % targets[0])
                inv[targets[0]] = LinComb({k: v.coeff(targets[0])}) # should invert coefficient, but mod 3 it doesn't matter
            else: # k --> a sum of multiple terms. Hope all but one term is already in the inverse mapping so far.
                old = [] # terms in orig[k] we already know how to write as orig[something]
                new = [] # terms in orig[k] we haven't seen before
                for val in targets:
                    if exclude_fn != None and exclude_fn(val):
                        continue
                    if val in inv.keys():
                        old.append(val)
                    else:
                        new.append(val)
                if len(new) == 0:
                    continue
                # k -> coeff * new[0] + sum coeff_i * old_i for old_i in old   implies
                # 1/coeff * (k - sum coeff_i * inv(old_i)) -> k
                elif len(new) == 1:
                    new_val = new[0]
                    new_keys = LinComb({k:1})
                    for val in old:
                        old_keys = inv[val]
                        new_keys.add_inplace(old_keys.scalar(-v.coeff(val)))
                    new_keys = new_keys.scalar(v.coeff(new[0])) # should be 1/v.coeff but same mod 3
                    inv[new_val] = new_keys
                else:
                    print("stuck on %s" % k)
                    for val in new:
                        inv[val] = LinComb({k:"?"})
        return inv


    def parse_mult(self, filename):
        rtn = {}
        fh = open(filename, "r")
        for line in fh:
            src, target = re.split("	->	", line)
            target = re.sub("\n", "", target)
            target_lc = LinComb.zero()
            for t in re.split("\+", target):
                if t != "o":
                    target_lc.add_inplace(LinComb({t:1}))
            rtn[src] = target_lc
        fh.close()
        return rtn

    def make_a0table(self):
        return self.parse_mult(self.file_prefix + "_BPAANSS_a0.txt")

    def make_h0table(self):
        return self.parse_mult(self.file_prefix + "_BPAANSS_h0.txt")

    def make_boc_h0table(self):
        return self.parse_mult(self.file_prefix + "_BPBocSS_h0.txt")

    def make_B2A(self):
        return self.parse_mult(self.file_prefix + "_BPB2A_table.txt")

    def alpha1(self, elt):
        # Input can be a string or LinComb; output is a LinComb
        if type(elt) == str:
            return self.h0table[elt]
        return elt.map(lambda s : self.h0table[s])

    def three(self, elt):
        # Input can be a string or LinComb; output is a LinComb
        if type(elt) == str:
            return self.a0table[elt]
        return elt.map(lambda s : self.a0table[s])

    def delta(self, elt):
        # Input is a string
        if not elt in self.bockstein.keys():
            raise Exception("%s not found as a Bockstein element" % elt)
        if self.bockstein[elt] == "":
            return LinComb.zero()
        rtn = self.bockstein[elt] # this is the actual answer; the rest is just name translation
        if rtn == "": # bockstein outputs strings, not LinComb's
            rtn = LinComb.zero()
        else:
            rtn = LinComb({rtn:1})

        def is_a0_divisible(lc):
            # a linear combination is 3-divisible if every generator is
            for t in lc.elts.keys():
                if not t in self.boc_a0div.keys():
                    return False
            return True

        # Can only apply B2A to non-3-divisible names, so we remove the
        # 3-multiples one-at-a-time, counting as we go.
        a0_count = 0
        while is_a0_divisible(rtn):
            rtn = rtn.map(lambda s : self.boc_a0div[s])
            a0_count += 1
        # now we can convert using B2A
        base = rtn.map(lambda s : self.B2A[s]) # answer / 3^{r-1}
        for i in range(a0_count - 1): # now put all those 3-multiples back in
            base = base.map(lambda s : self.boc_a0table[s])
        return base

    def beta_n(self, n, elt):
        # Input is a string or LinComb
        # n = "1" ==> multiply by beta1 (determined by self.beta_mapping)
        if type(elt) == str:
            rtn = self.B2A_inv[elt]
            rtn = rtn.map(lambda s : self.beta_table['beta' + n][s])
            rtn = rtn.map(lambda s : self.delta(s))
            return rtn
        return elt.map(lambda s : self.beta_n(n, s))

    def beta1(self, elt, n=1):
        if n == 1:
            return self.beta_n("1", elt)
        else:
            return self.beta1(self.beta1(elt),n-1)




    def beta2(self, elt):
        return self.beta_n("2", elt)
    def beta33(self, elt):
        return self.beta_n("33", elt)
    def beta4(self, elt):
        return self.beta_n("4", elt)
    def beta5(self, elt):
        return self.beta_n("5", elt)
    def beta63(self, elt):
        return self.beta_n("63", elt)

    def in_deg(self, x, y):
        return [i for i,r in self.deg_table.items() if r['s'] == x and r['f'] == y]

    def show_deg(self, name):
        return (self.deg_table[name]['s'], self.deg_table[name]['f'], self.deg_table[name]['nov'])

