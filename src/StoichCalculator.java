import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

// Use apache math library for fraction and lowest common 
// multiple classes and methods to make coefficients integers
import org.apache.commons.math3.fraction.BigFraction;
import org.apache.commons.math3.util.ArithmeticUtils;

public class StoichCalculator {

	// Attempts to balance equation but has three possible scenarios: one in
	// which there is an error in balancing because
	// the equation is impossible to balance which is represented by a return
	// value of 1; one in which the equation is already balanced
	// correctly which is represented by a return value of 0; and one in which
	// the equation is possible to balance but incorrectly balanced
	// which is corrected by the program. I use linear algebra to solve for the
	// balanced equation using rref, methods have been referenced from
	// a combination of online resources because I was not entirely familiar
	// with linear algebra.

	public static int balanceEquation(ArrayList<Compound> reactants, ArrayList<Compound> products,
			HashMap<Integer, Integer> listOfElements) {

		// Create list of all elements
		ArrayList<Compound> equationList = new ArrayList<Compound>();
		equationList.addAll(reactants);
		equationList.addAll(products);

		double matrix[][] = new double[listOfElements.size()][equationList.size()]; // check
		double coefficientsDouble[] = new double[equationList.size()];
		int coefficientsFinal[] = new int[equationList.size()];

		// Create equations for each element based on coefficients and store
		// them in matrix
		for (int i = 0; i < equationList.size(); i++) {
			ArrayList<Element> e = equationList.get(i).getElementList();
			for (int j = 0; j < e.size(); j++) {
				int id = e.get(j).getId();
				matrix[listOfElements.get(id)][0] += e.get(j).getSubscript() * equationList.get(i).getCoefficient();
			}
		}

		// Refer above for return values explanation, determines if equation is
		// possible to balance but incorrectly balanced, if it is already
		// balanced or if it is impossible to balance
		if (isBalance(matrix, reactants.size())) {
			return 0;
		} else {
			// Use rref to solve for the cofficients of the equation
			matrix = determineReducedRowEchelonForm(matrix);

			// Extract the coefficients from matrix
			assignCoefficients(coefficientsDouble, matrix);

			// Make coefficients integers since you cannot have half a molecule
			coefficientsFinal = finalizeCoefficients(coefficientsDouble);

			// Set coefficients of equation to the determined coefficients

			for (int i = 0; i < reactants.size(); i++) {
				reactants.get(i).setCoefficient(coefficientsFinal[i]);
			}

			for (int i = 0; i < products.size(); i++) {
				products.get(i).setCoefficient(coefficientsFinal[i + reactants.size()]);
			}

			// Checks if matrix is impossible to balance by checking if is
			// balanced after being balanced, returns appropriate values (see
			// above)
			matrix = new double[listOfElements.size()][equationList.size()];
			for (int i = 0; i < equationList.size(); i++) {
				ArrayList<Element> e = equationList.get(i).getElementList();
				for (int j = 0; j < e.size(); j++) {
					int id = e.get(j).getId();
					matrix[listOfElements.get(id)][i] += e.get(j).getSubscript() * equationList.get(i).getCoefficient();
				}
			}
			if (!isBalance(matrix, reactants.size()))
				return 1;
		}
		return 2;
	}

	// Uses matrix of elements in equation to check if it is balanced, checks if
	// there are the same number of that element on the reactants side as there
	// are on the products side
	public static boolean isBalance(double matrix[][], int buffer) {
		double r[] = new double[matrix.length], p[] = new double[matrix.length];

		for (int i = 0; i < buffer; i++) {
			for (int j = 0; j < r.length; j++) {
				r[j] += matrix[j][i];
			}
		}

		for (int i = buffer; i < matrix[0].length; i++) {
			for (int j = 0; j < p.length; j++) {
				p[j] += matrix[j][i];
			}
		}

		for (int i = 0; i < r.length; i++) {
			if (r[i] != p[i])
				return false;
		}
		return true;
	}

	// Assigns coefficient from the matrix to the coefficient array because you
	// need to extract them from row reduced format
	public static void assignCoefficients(double[] coefficients, double matrix[][]) {

		for (int i = 0; i < matrix.length; i++) {
			coefficients[i] = Math.abs(matrix[i][matrix[i].length - 1]);
		}

		for (int i = matrix.length; i < coefficients.length; i++) {
			coefficients[i] = 1;
		}
	}

	// Returns the integer format of double coefficients, by taking the fraction
	// form and multiplying by the lowest common multiple of the denominators
	public static int[] finalizeCoefficients(double[] coefficients) {
		BigFraction rationalCoefficients[] = new BigFraction[coefficients.length];

		// Makes each double into fraction
		for (int i = 0; i < coefficients.length; i++) {
			rationalCoefficients[i] = new BigFraction(coefficients[i], 0.00000002D, 10000);
			System.out.println(rationalCoefficients[i].toString() + "  ");
		}

		// Gets lowest common multiple of denominators
		int lcm = Integer.MIN_VALUE;
		for (int i = 0; i < coefficients.length - 1; i++) {
			lcm = Math.max(ArithmeticUtils.lcm(rationalCoefficients[i].getDenominatorAsInt(),
					rationalCoefficients[i + 1].getDenominatorAsInt()), lcm);
		}

		// Finally converts all fractions into integers by multiplying by lcm of
		// denominators
		int coefficientsFinal[] = new int[coefficients.length];
		for (int i = 0; i < coefficients.length; i++) {
			coefficientsFinal[i] = (int) (coefficients[i] * lcm);
		}

		return coefficientsFinal;
	}

	// Determines quantity of all elements possible using available data, or if
	// there is insufficient data
	public static boolean determineQuantity(ArrayList<Compound> reactants, ArrayList<Compound> products) {
		ArrayList<Compound> equationList = new ArrayList<Compound>();
		equationList.addAll(reactants);
		equationList.addAll(products);

		boolean insufficientData = true;

		// Goes through pecking order of required information to be able
		// determine ratios and mass for all compounds
		for (int i = 0; i < equationList.size(); i++) {
			Compound c = equationList.get(i);
			if (c.getMoles() > 0) {
				c.setMass(c.getMoles() * c.getMolarMass());
				insufficientData = false;
			} else if (c.getMass() > 0) {
				c.setMoles(c.getMass() / c.getMolarMass());
				insufficientData = false;
			} else if (c.getVolume() > 0 && c.getDensity() > 0) {
				c.setMass((c.getVolume() * c.getDensity()) / c.getMolarMass());
				c.setMoles(c.getMass() / c.getMolarMass());
				insufficientData = false;
			}
		}

		// If sufficient data not determined return false
		if (insufficientData)
			return false;
		else
			return true;
	}

	// Determine the limiting reagent of the equation
	public static int determineLimitingReagent(ArrayList<Compound> reactants) {
		int index = 0;
		double min = Double.MAX_VALUE;

		// Determines index of limiting reagent by going through the reactants
		// and determing the one with the smallest quantity / coefficient
		for (int i = 0; i < reactants.size(); i++) {
			if (reactants.get(i).getMoles() / reactants.get(i).getCoefficient() < min
					&& reactants.get(i).getMoles() != 0) {
				min = reactants.get(i).getMoles() / reactants.get(i).getCoefficient();
				index = i;
			}
		}
		return index;
	}

	// Determines all the unknowns of the reactants and products based on
	// limiting reagent
	public static void determineUnknowns(ArrayList<Compound> reactants, ArrayList<Compound> products,
			int limitingReagentIndex) {

		Compound limitingReagent = reactants.get(limitingReagentIndex);

		// Calculates the moles and mass of all reactants by mutltipying the
		// limiting reagent by the stoichiometric ratios between the two
		for (int i = 0; i < reactants.size(); i++) {
			Compound c = reactants.get(i);
			if (c.getMoles() == 0 && c.getMass() == 0) {
				c.setMoles(
						((double) c.getCoefficient() / limitingReagent.getCoefficient()) * limitingReagent.getMoles());
				c.setMass(c.getMoles() * c.getMolarMass());
			}
		}

		// Calculates the moles and mass of all products by mutltipying the
		// limiting reagent by the stoichiometric ratios between the two
		for (int i = 0; i < products.size(); i++) {
			Compound c = products.get(i);
			c.setMoles(((double) c.getCoefficient() / limitingReagent.getCoefficient()) * limitingReagent.getMoles());
			c.setMass(c.getMoles() * c.getMolarMass());
		}

	}

	//The code below is from other sources, it is the code for the linear algebra to determine the balanced equation
	
	public static double[][] determineReducedRowEchelonForm(double[][] matrix) {
		int lead = 0;
		int rowCount = matrix.length;
		int colCount = matrix[0].length;
		int i;
		boolean quit = false;

		for (int row = 0; row < rowCount && !quit; row++) {
			if (colCount <= lead) {
				quit = true;
				break;
			}

			i = row;

			while (!quit && matrix[i][lead] == 0) {
				i++;
				if (rowCount == i) {
					i = row;
					lead++;

					if (colCount == lead) {
						quit = true;
						break;
					}
				}
			}

			if (!quit) {
				swapRows(matrix, i, row);

				if (matrix[row][lead] != 0)
					multiplyRow(matrix, row, 1.0f / matrix[row][lead]);

				for (i = 0; i < rowCount; i++) {
					if (i != row)
						subtractRows(matrix, matrix[i][lead], row, i);
				}
			}
		}
		return matrix;
	}

	public static void swapRows(double[][] matrix, int rowIndex1, int rowIndex2) {
		double[] swap = new double[matrix[0].length];

		for (int c1 = 0; c1 < matrix[0].length; c1++)
			swap[c1] = matrix[rowIndex1][c1];

		for (int c1 = 0; c1 < matrix[0].length; c1++) {
			matrix[rowIndex1][c1] = matrix[rowIndex2][c1];
			matrix[rowIndex2][c1] = swap[c1];
		}
	}

	public static void multiplyRow(double[][] matrix, int rowIndex, double scalar) {
		for (int c1 = 0; c1 < matrix[0].length; c1++)
			matrix[rowIndex][c1] *= scalar;
	}

	public static void subtractRows(double[][] matrix, double scalar, int subtractScalarFromRowIndex, int rowIndex) {
		for (int c1 = 0; c1 < matrix[0].length; c1++)
			matrix[rowIndex][c1] -= scalar * matrix[subtractScalarFromRowIndex][c1];
	}
}
