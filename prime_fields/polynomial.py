"""
Define the Polynomial class that can be applied to form companion matrices.
"""
from numbers import Number, Integral
from primes import Matrix, Element


class Polynomial:
    """Instantiate a polynomial with coefficients and degree."""
    def __init__(self, coefs):
        if len(coefs) == 0:
            self.coefficients = (0, )
        else:
            coefi = list(coefs)
            while coefi[-1] == 0 and len(coefi) != 1:
                coefi.pop()
            self.coefficients = tuple(coefi)

    @property
    def degree(self):
        return len(self.coefficients) - 1

    def __str__(self):
        coefs = self.coefficients
        terms = []
        if coefs[0]:
            terms.append(str(coefs[0]))
        if self.degree and coefs[1]:
            terms.append(f"{'' if coefs[1] == 1 else coefs[1]}x")
        terms += [f"{'' if c == 1 else c}x^{d}"
                  for d, c in enumerate(coefs[2:], start=2) if c]
        return " + ".join(reversed(terms)) or "0"

    def __repr__(self):
        return self.__class__.__name__ + "(" + repr(self.coefficients) + ")"

    def __eq__(self, other):
        return isinstance(other, Polynomial) and\
             self.coefficients == other.coefficients

    def __add__(self, other):
        if isinstance(other, Polynomial):
            common = min(self.degree(), other.degree()) + 1
            coefs = tuple(a + b for a, b in zip(self.coefficients,
                                                other.coefficients))
            coefs += self.coefficients[common:] + other.coefficients[common:]
            return Polynomial(coefs)
        elif isinstance(other, Number):
            return Polynomial((self.coefficients[0] + other,)
                              + self.coefficients[1:])
        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Polynomial):
            maxi = max(self.degree(), other.degree())+1
            coefs = []
            for i in range(maxi):
                if i < len(self.coefficients):
                    coefs_self = self.coefficients[i]
                else:
                    coefs_self = 0
                if i < len(other.coefficients):
                    coefs_other = other.coefficients[i]
                else:
                    coefs_other = 0
                coefs.append(coefs_self - coefs_other)
            return Polynomial(tuple(coefs))
        elif isinstance(other, Number):
            return Polynomial((self.coefficients[0] - other,)
                              + self.coefficients[1:])
        else:
            raise NotImplementedError

    def __rsub__(self, other):
        if isinstance(other, Polynomial):
            maxi = max(self.degree(), other.degree())+1
            coefs = []
            for i in range(maxi):
                if i < len(self.coefficients):
                    coefs_self = self.coefficients[i]
                else:
                    coefs_self = 0
                if i < len(other.coefficients):
                    coefs_other = other.coefficients[i]
                else:
                    coefs_other = 0
                coefs.append(coefs_other - coefs_self)
            return Polynomial(tuple(coefs))
        elif isinstance(other, Number):
            coefs = [other - self.coefficients[0]]
            for i in range(1, len(self.coefficients)):
                coefs.append(-self.coefficients[i])
            return Polynomial(tuple(coefs))
        else:
            raise NotImplementedError

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            deg = self.degree() + other.degree()
            coefs = [0 for i in range(deg + 1)]
            for i in range(self.degree() + 1):
                for j in range(other.degree() + 1):
                    coefs[i+j] += self.coefficients[i] * other.coefficients[j]
            return Polynomial(tuple(coefs))
        elif isinstance(other, Number):
            coefs = [other * i for i in self.coefficients]
            return Polynomial(tuple(coefs))
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self * other

    def __pow__(self, other):
        if isinstance(other, Integral):
            poly, ite = self, self
            for i in range(1, other):
                poly = poly * ite
            return poly
        else:
            raise NotImplementedError

    def __call__(self, other):
        if isinstance(other, Number):
            val = 0
            for i in range(self.degree() + 1):
                val += self.coefficients[i] * (other**i)
            return val
        else:
            raise NotImplementedError

    def dx(self):
        if self.degree() >= 1:
            return Polynomial(tuple(self.coefficients[i]*i
                                    for i in range(1, self.degree() + 1)))
        else:
            return Polynomial((0,))


def derivative(poly):
    return poly.dx()


def companion_matrix(poly):
    if not isinstance(poly, Polynomial):
        raise TypeError("Input needs to be a Polynomial"
                        "with at least degree 1.")
    elif not poly.degree - 1:
        raise ValueError("Polynomial needs to have a degree of at least 1.")
    n = poly.degree
    coefs = [coefficient/poly.coefficients[-1]
             for coefficient in poly.coefficients]
    coefs.pop()
    if isinstance(poly.coefficients[0], Element):
        field = poly.coefficients[0][0].field
    else:
        field = 0
    companion = []
    for i in range(n):
        row_vector = [0 for j in range(n)]
        if i:
            row_vector[i - 1] = 1
        row_vector[-1] = -coefs.pop()
        companion.append(row_vector)
    return Matrix(n, companion, field)


def block_diagonal(square_matrices):
    """Create a block_diagonal matrix with the given square matrices."""
    field = _check_input(square_matrices)
    total_size = sum(i.size for i in square_matrices)
    current_size = 0
    blockdiag_matrix = []
    for matrix in square_matrices:
        blockdiag_matrix.extend(_populate_zeroes(matrix,
                                                 current_size, total_size))
        current_size += matrix.size
    return Matrix(total_size, blockdiag_matrix, field)


def _check_input(square_matrices):
    for matrix in square_matrices:
        if not isinstance(matrix, Matrix):
            raise TypeError("Not all items in list are matrices.")
        m = square_matrices[0].field
        if not matrix.field == m:
            raise TypeError("Not all matrices are defined on same field.")
    return square_matrices[0].field


def _populate_zeroes(matrix, current_size, total_size):
    rows = matrix.size
    resulting_rows = [[0 for i in range(current_size)] for row in range(rows)]
    for i in range(rows):
        resulting_rows[i].extend(matrix.matrix[i])
        resulting_rows[i].extend([0 for j in
                                  range(total_size - current_size - rows)])
    return resulting_rows
