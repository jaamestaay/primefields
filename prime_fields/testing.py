import polynomial as poly
import primes

F2 = primes.Fp(2)
a_1 = poly.companion_matrix(poly.Polynomial([1, 1]), 2)
print(a_1**3 + poly.identity(1, F2))
a_2 = poly.companion_matrix(poly.Polynomial([1, 1])**2, 2)
print(a_2**3 + poly.identity(2, F2))
a_3 = poly.companion_matrix(poly.Polynomial([1, 1, 1])**3, 2)
print(a_3**3 + poly.identity(6, F2))
a_4 = poly.companion_matrix(poly.Polynomial([1, 1, 0, 0, 1])**2, 2)
print(a_4**3 + poly.identity(8, F2))
a = poly.block_diagonal([a_1, a_2, a_3, a_4])
print(a)
expr_a = a**3 + poly.identity(17, F2)
print(expr_a)
