import copy
import frozendict

class LinComb(object):
    # elts is a dict whose keys are strings not containing the null character
    # and values are ints (interpreted as in F_3) or "?"
    prime = 3

    def __init__(self, elts={}):
        elts = { k:v for k,v in elts.items() if v == "?" or v % self.__class__.prime != 0 }
        self.elts = elts

    def __repr__(self):
        if self.elts == {}:
            return "0"
        rtn = ""
        for elt,coeff in self.elts.items():
            if coeff == "?":
                coeff_str = "?"
            elif coeff % self.__class__.prime == 1:
                coeff_str = ""
            elif coeff % self.__class__.prime != 0:
                coeff_str = str(coeff)
            else:
                raise Exception("coefficient of %s is zero" % elt)
            # assume coeff != 0
            if rtn == "":
                rtn += "%s%s" % (coeff_str, elt)
            else:
                rtn += " + %s%s" % (coeff_str, elt)
        return rtn

    # Needed in order to make dicts with LinComb objects as keys
    def __hash__(self):
        return hash(frozendict.frozendict(self.elts))

    def serialize(self):
        rtn = ""
        for k,coeff in self.elts.items():
            rtn += "%s\1%s\1" % (coeff, k)
        return rtn

    @classmethod
    def deserialize(cls, s):
        if s == "":
            return cls({})
        elts = {}
        lst = re.split("\1", s)
        # lst contains alternating coeff, generator pairs, with an "" at the
        # end since s is terminated by the null character
        for i in range(int((len(lst)-1)/2)):
            if lst[2*i] == "?":
                elts[lst[2*i+1]] = "?"
            else:
                elts[lst[2*i+1]] = int(lst[2*i])
        return cls(elts)


    def __bool__(self):
        return not self.iszero()

    def __eq__(self, other):
        return self.elts == other.elts

    def __ne__(self, obj):
            return not self == obj

    def iszero(self):
        return self.elts == {}

    def firstkey(self):
        # returns one
        if self.elts == {}:
            return None
        return list(self.elts.keys())[0]

    @classmethod
    def zero(cls):
        return cls({})

    def coeff(self, key):
        if not key in self.elts.keys():
            return 0
        else:
            return self.elts[key]

    def add(self, other):
        # returns self + other (without changing self)
        rtn = copy.deepcopy(self)
        for k,v in other.elts.items():
            if k in self.elts.keys():
                if v == "?" or self.elts[k] == "?":
                    rtn.elts[k] = "?"
                    continue
                coeff = (v + self.elts[k]) % self.__class__.prime
                if coeff == 0:
                    rtn.elts.pop(k)
                else:
                    rtn.elts[k] = coeff
            else:
                rtn.elts[k] = v
        return rtn

    def add_inplace(self, other):
        # sets self = self + other
        rtn = {}
        for k,v in other.elts.items():
            if k in self.elts.keys():
                if v == "?" or self.elts[k] == "?":
                    self.elts[k] = "?"
                    continue
                coeff = (v + self.elts[k]) % self.__class__.prime
                if coeff == 0:
                    self.elts.pop(k)
                else:
                    self.elts[k] = coeff
            else:
                self.elts[k] = v

    def scalar(self, c):
        # For int c (or c = ?), returns c * self (without changing self)
        rtn = copy.deepcopy(self)
        if c == "?":
            for k in self.elts.keys():
                rtn.elts[k] = "?"
        elif c % self.__class__.prime == 0:
            rtn.elts = {}
        else:
            for k,v in self.elts.items():
                # if v == "?" then it stays as "?"
                if v != "?":
                    rtn.elts[k] = (c * v) % self.__class__.prime
        return rtn

    def map(self, fn):
        # This applies fn: keys -> LinComb to each key in self, linearly.
        # Returns a new object without changing self.
        rtn = self.__class__({})
        for k,coeff in self.elts.items():
            newterm = fn(k)
            rtn.add_inplace(newterm.scalar(coeff))
        return rtn


