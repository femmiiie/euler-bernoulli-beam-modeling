import matplotlib.pyplot as plt
from board import Board
from functions import *


def step2():
	inorm = lambda i : max(abs(min(i)), max(i))

	n = 10
	E = 1.3 * 10**10
	free_board = Board(unweighted, y_correct, 9.81, 30, 3, 0, 2, E)
	free_board.set_n(n)

	x_vals = free_board.get_x_vals()

	x = free_board.solve()
	fe = free_board.get_fe()
	correct = free_board.get_correct()

	print(x_vals)
	print(correct)
	print(x)


	print(f"Forward error of step 1 calculation: {inorm(fe)}")

	plt.style.use('_mpl-gallery')
	fig, ax = plt.subplots()
	ax.scatter(x_vals, x, [15] * n, [[1, 0, 0]])
	ax.scatter(x_vals, fe, [15] * n, [[0, 1, 0]])
	ax.scatter(x_vals, correct, [15] * n, [[0, 0, 1]])

	plt.show()


def step3():
	E = 1.3 * 10**10
	free_board = Board(unweighted, y_correct, 9.81, 30, 3, 0, 2, E)
	free_board.pow_2_chart(11)


def step5():
	E = 1.3 * 10**10
	sinu_board = Board(sinusoidal, y_correct, 9.81, 30, 3, 100, 2, E)
	sinu_board.pow_2_chart(11)


def step5_5():
	E = 1.3 * 10**10
	L = 2
	max_k = 11

	h_vals = []
	error_vals = []

	for k in range(1, max_k + 1):
		n = 10 * 2**k
		h = L / n
		h_vals.append(h)

		board = Board(sinusoidal, y_correct, 9.81, 30, 3, 100, L, E)
		board.set_n(n)

		x_numerical = board.solve()
		x_vals = board.get_x_vals()

		x_L = x_vals[-1]
		y_L_numerical = x_numerical[-1]
		y_L_exact = y_correct(x_L, board.f(x_L), 30, 3, L, E, board.I)

		error = abs(y_L_numerical - y_L_exact)
		error_vals.append(error)

	fig, ax = plt.subplots()
	ax.set_xscale("log", base=10)
	ax.set_yscale("log", base=10)

	plt.figure(figsize=(8, 6))
	plt.loglog(h_vals, error_vals, 'o-')
	plt.xlabel('h (grid spacing)')
	plt.ylabel('Absolute Error at x=L')
	plt.title('Error vs h on log-log scale (Sinusoidal Load)')
	plt.grid(True, which="both", ls="--")
	plt.show()

	print(h_vals)
	print(error_vals)


def step6():
	E = 1.3 * 10**10
	diver_board = Board(diver, y_correct, 9.81, 30, 3, 100, 2, E)
	diver_board.set_n(1280) #optimal value from step 5
	vals = diver_board.solve()

	for val in vals:
		print(val)


step2()
step3()
step5()
step5_5()
step6()
