import java.util.ArrayList;
import java.util.HashMap;

public class ChemicalEquation {
	private ArrayList<Compound> reactants = new ArrayList<Compound>();
	private ArrayList<Compound> products = new ArrayList<Compound>();

	// Hashmap of each element in equation, Hashmap is necessary for matrix in balancing equation
	private HashMap<Integer, Integer> listOfElements;

	// boolean which checks if it is technically possibly to balance the equation
	private boolean canBeBalanced = true;

	// Constructor initializes instance variables
	public ChemicalEquation(ArrayList<Compound> r, ArrayList<Compound> p) {
		this.reactants = r;
		this.products = p;
		this.listOfElements = createUniqueElementList(r, p);
	}

	// Initializes hashmap, creating the list of elements (no duplicates)
	public HashMap<Integer, Integer> createUniqueElementList(ArrayList<Compound> r, ArrayList<Compound> p) {
		HashMap<Integer, Integer> listOfElements = new HashMap<Integer, Integer>();

		// Boolean to check if element has been added to the hashmap
		boolean visited[] = new boolean[119];

		int counter = 0;

		// Makes a hashmap of all unique elements of the chemical equation and returns it to instance variable
		for (int i = 0; i < r.size(); i++) {
			Compound compound = r.get(i);

			// Try/Catch for a specific impossible to balance equation such as
			// incomplete combustion is checked, this makes sense programmatically
			try {
				for (int j = 0; j < compound.getElementList().size(); j++) {
					Element element = compound.getElementList().get(j);
					if (!visited[element.getId()]) {
						listOfElements.put(element.getId(), counter++);
						visited[element.getId()] = true;
					}
				}
			} catch (NullPointerException e) {
				canBeBalanced = false;
				return null;
			} catch (ArrayIndexOutOfBoundsException e) {
				canBeBalanced = false;
				return null;
			}

		}
		return listOfElements;
	}

	public boolean solveEquation() {
		
		// If determined to be impossible to balance in createUniqueElementList method or checkElementQuantityInput methods, print error 
		if (!canBeBalanced) {
			Console.messageList.add("Cannot balance, Retry");
			return false;
		} else if (!checkElementQuantityInput(this.reactants, this.products)) {
			Console.messageList.add("I can't balance this! Try Again!");
			return false;
		}

		// Attempts to balance equation, calls calculation method, refer to method for specifics
		int balanceCheck = StoichCalculator.balanceEquation(this.reactants, this.products, this.listOfElements);

		// If equation determined to be impossible to balance in process, print error
		// Else if equation is incorrectly balanced, print that fact but still balance it
		if (balanceCheck == 1) {
			Console.messageList.add("I can't balance this! Try Again!");
			return false;
		} else if (balanceCheck == 2) {
			EquationPanel.compoundList.clear();
			EquationPanel.compoundList.addAll(this.reactants);
			EquationPanel.compoundList.addAll(this.products);
			Console.messageList.add("Incorrectly balanced! I had to do it for you!");
		}

		// Runs stoichiometry methods to determine the ratios, will return error if insufficient data  
		if (StoichCalculator.determineQuantity(this.reactants, this.products)) {
			// Gets limiting reagent
			int limitingReagentIndex = StoichCalculator.determineLimitingReagent(this.reactants);
			
			// Determines unknowns using limiting reagent
			StoichCalculator.determineUnknowns(this.reactants, this.products, limitingReagentIndex);
			EquationPanel.compoundList.clear();
			EquationPanel.compoundList.addAll(this.reactants);
			EquationPanel.compoundList.addAll(this.products);
		} else {
			Console.messageList.add("You didn't give me enough data! Try Again!");
			return false;
		}

		return true;
	}

	// Checks for potential error inputs that can be easily detected such as empty compounds or mistmatched input
	public boolean checkElementQuantityInput(ArrayList<Compound> reactants, ArrayList<Compound> products) {
		
		ArrayList<Integer> reactantElements = new ArrayList<Integer>();
		ArrayList<Integer> productElements = new ArrayList<Integer>();

		// Adds all individual elements id's from reactants into a list
		for (int i = 0; i < reactants.size(); i++) {
			// Checks if compound is empty or if coefficient is 0 and returns error since that is invalid
			if (reactants.get(i).getElementList().size() == 0 || reactants.get(i).getCoefficient() == 0) {
				return false;
			}

			for (int j = 0; j < reactants.get(i).getElementList().size(); j++) {
				reactantElements.add(reactants.get(i).getElementList().get(j).getId());
			}
		}

		// Adds all individual elements id's from products into a list
		for (int i = 0; i < products.size(); i++) {
			// Checks if compound is empty or if coefficient is 0 and returns error since that is invalid
			if (products.get(i).getElementList().size() == 0 || products.get(i).getCoefficient() == 0) {
				return false;
			}
			
			for (int j = 0; j < products.get(i).getElementList().size(); j++) {
				productElements.add(products.get(i).getElementList().get(j).getId());
			}
		}

		// Makes sure all elements in reactants appear in products
		for (int i = 0; i < reactantElements.size(); i++) {
			if (!productElements.contains(reactantElements.get(i))) {
				return false;
			}
		}

		// Makes sure all elements in products appear in reactants
		for (int i = 0; i < productElements.size(); i++) {
			if (!reactantElements.contains(productElements.get(i))) {
				return false;
			}
		}
		return true;
	}
}
