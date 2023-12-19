"""
Define the finite field of prime p elements Fp, to be used in matrix algebra.
"""
from functools import wraps
from numbers import Integral
from sympy import isprime


class Fp:
    """Define the Field for which the matrices are defined on."""
    def __init__(self, prime):
        if isprime(prime):
            self.prime = prime
        else:
            raise ValueError("Value is not prime.")

    def __str__(self):
        return f"F{self.prime}"

    def __repr__(self):
        return f"{type(self).__name__}({self.prime!r})"

    def __eq__(self, other):
        if isinstance(other, Integral):
            return self.prime == other
        elif isinstance(other, Fp):
            return self.prime == other.prime
        else:
            raise NotImplementedError("Equality has not been implemented.")


def make_other_Element(meth):
    """Check if other is a numbers.Integral and change value accordingly."""
    @wraps(meth)
    def fn(self, other):
        if isinstance(other, Element):
            if self.field.prime != other.field.prime:
                raise TypeError("Other must be an Element of same field.")
        elif isinstance(other, Integral):
            other = Element(self.field, other)
        else:
            raise TypeError("Other must be an integer or an Element.")
        return meth(self, other)
    return fn


class Element:
    """Define the Element within a matrix in the field Fp."""
    def __init__(self, field, value):
        if not isinstance(field, Fp):
            raise TypeError("Define field using Fp.")
        self.field = field
        if not isinstance(value, Integral):
            raise TypeError("Value should be an integer.")
        self.value = value % field.prime

    def __str__(self):
        return f"{self.value} mod {self.field.prime}"

    # def __repr__(self):
    #   return f"{type(self).__name__}({self.field}, {self.value})"
    # In the interest of making the matrix readable, I've opted to define
    # __repr__ to give just self.value. The "correct" implementation of
    # __repr__ is made as notes above.
    def __repr__(self):
        return f"{self.value}"

    @make_other_Element
    def __add__(self, other):
        return Element(self.field, self.value + other.value)

    def __radd__(self, other):
        return self + other

    @make_other_Element
    def __sub__(self, other):
        return Element(self.field, self.value - other.value)

    @make_other_Element
    def __rsub__(self, other):
        return Element(self.field, other.value - self.value)

    @make_other_Element
    def __mul__(self, other):
        return Element(self.field, self.value * other.value)

    def __rmul__(self, other):
        return self * other

    # __pow__ will only be defined for integer values, not Element.
    def __pow__(self, other):
        if isinstance(other, Integral):
            return Element(self.field, self.value ** other)
        else:
            raise TypeError("Exponential should be an integer.")


def _validate(array):
    size = len(array)
    for i in range(size):
        if not len(array[i]) is size:
            return False
    return True


def _elementarise(field, array):
    if not field:
        return array
    n = len(array)
    for i in range(n):
        for j in range(n):
            if isinstance(array[i][j], Element):
                if array[i][j].field.prime == field:
                    continue
                else:
                    # Slightly dangerous implementation, might change value.
                    # But just covers the fool case.
                    array[i][j] = Element(field, array[i][j].value)
            elif not isinstance(array[i][j], Integral):
                raise TypeError("All matrix elements should be integers.")
            else:
                array[i][j] = Element(field, array[i][j])
    return array


class Matrix:
    """Define the square matrix with Element elements in the field Fp."""
    def __init__(self, array, field=0):
        if not field:
            pass
        elif not isinstance(field, Fp):
            raise TypeError("Define field using Fp.")
        elif not _validate(array):
            raise TypeError("Array does not give a square matrix.")
        self.size = len(array)
        self.field = field
        self.matrix = _elementarise(field, array)

    def __str__(self):
        d = "\n ".join(str(submatrix) for submatrix in self.matrix).__str__()
        return f"[{d}]"

    # def __repr__(self):
    #    return f"{type(self).__name__}"\
    #           f"({self.field}, {self.size},{self.matrix})"
    # Yet again, by design choice, opting for a similar __repr__ as __str__.
    # "Correct" implementation is made as notes above.
    def __repr__(self):
        d = "\n ".join(str(submatrix) for submatrix in self.matrix).__str__()
        return f"[{d}]"

    def __add__(self, other):
        """Add two matrices together element-wise."""
        _checkmatrix(self, other)
        result_matrix = []
        for i in range(self.size):
            result_matrix.append(list())
            for j in range(self.size):
                result_matrix[i].append(self.matrix[i][j] + other.matrix[i][j])
        return Matrix(result_matrix, self.field)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        _checkmatrix(self, other)
        result_matrix = []
        for i in range(self.size):
            result_matrix.append(list())
            for j in range(self.size):
                result_matrix[i].append(self.matrix[i][j] - other.matrix[i][j])
        return Matrix(result_matrix, self.field)

    def __mul__(self, other):
        result_matrix = []
        if isinstance(other, Integral):
            for i in range(self.size):
                result_matrix.append(list())
                for j in range(self.size):
                    result_matrix[i].append(self.matrix[i][j] * other)
        elif isinstance(other, Matrix):
            for i in range(self.size):
                result_matrix.append(list())
                for j in range(self.size):
                    new_value = sum(self.matrix[i][k] * other.matrix[k][j]
                                    for k in range(self.size))
                    result_matrix[i].append(new_value)
        return Matrix(result_matrix, self.field)

    def __rmul__(self, other):
        return self * other

    def __pow__(self, other):
        if not isinstance(other, Integral):
            raise TypeError("Exponent must be an integer.")
        result_matrix = self
        for i in range(other - 1):
            result_matrix *= self
        return result_matrix


def _checkmatrix(self, other):
    if not isinstance(other, type(self)):
        raise TypeError("At least one of the operation" /
                        "elements is not a Matrix.")
    elif self.size - other.size:
        raise TypeError("Matrices must be the same size.")
    elif not self.field or self.field.prime - other.field.prime:
        raise TypeError("Matrices must be defined in the same field.")
