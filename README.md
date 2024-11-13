# BinPackingProblem

1. Instances

   	You can find instances in https://drive.google.com/drive/folders/1Wk_TFV60XMFpMbwzrVNddiPPpYjF4axb?usp=sharing

	The file names follow the pattern:

		First number: class
		Second number: quantity
		Last number: instance index

	Classes from 1 to 8 are the 8 classes of Martello, Pisinger e Vigo (2000) problems.
	Class 9 is Ivancic, Mathur e Mohanty (1989) problems
	Class 10 is Thapatsuwan et al. (2012) problems

3. Residual Space and Orientations rules quantity parameters

	Classes from 1 to 8
		Residual Space rules quantity 2
		Orientation rules quantity 1

	Classes from 9 to 10
		Residual Space rules quantity 7
		Orientation rules quantity 6

4. Layout results.csv file:
	Column 1: Classe
	Column 2: Quantity
	Column 3: Instance Index
	Column 4: Objective Function Value
	Column 5: time in second
	Column 6: waste (used only in class 10 problems)
	Column 7: solution: layout [[[idItem, xi, yi, zi, orientation],...] 
