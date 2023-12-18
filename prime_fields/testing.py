import polynomial as poly
import primes

F2 = primes.Fp(2)
a_1 = poly.companion_matrix(poly.Polynomial([1, 1]), 2)
a_2 = poly.companion_matrix(poly.Polynomial([1, 1])**2, 2)
a_3 = poly.companion_matrix(poly.Polynomial([1, 1, 1])**3, 2)
a_4 = poly.companion_matrix(poly.Polynomial([1, 1, 0, 0, 1])**2, 2)
a = poly.block_diagonal([a_1, a_2, a_3, a_4])
print(a)
expr_a = a**3 + poly.identity(17, F2)
print(expr_a)
