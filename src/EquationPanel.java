import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.JTextField;

public class EquationPanel extends JPanel implements MouseListener, KeyListener {

	private static final long serialVersionUID = 4631820944519889788L;

	// The index of the current compound that is selected (used for choosing
	// which data/formula panel to display)
	static int compoundIndex = 0;

	static Images background, addReactant, subReactant, addProduct, subProduct, plus;

	// Lists of reactant and product images
	static ArrayList<Images> reactantBox = new ArrayList<Images>();
	static ArrayList<Images> productBox = new ArrayList<Images>();

	// List of textboxes for the coefficients of the compounds
	static ArrayList<JTextField> coefficientList = new ArrayList<JTextField>();
	static ArrayList<Compound> compoundList = new ArrayList<Compound>();

	// Holds the current font size
	static int fontFactor = 7;

	// Calls the constructor, which creates a new blank equation
	EquationPanel() throws IOException {
		resetLists();
		this.addMouseListener(this);
		this.addKeyListener(this);
		this.requestFocus();

		background = new Images("EquationPanel.png", 0, 0);

		reactantBox.add(new Images("reactant.png", 312, 275));
		coefficientList.add(new JTextField(2));
		this.add(coefficientList.get(0));
		coefficientList.get(0).setBounds(reactantBox.get(0).x + reactantBox.get(0).img.getWidth() / 2 - 15,
				reactantBox.get(0).y + reactantBox.get(0).img.getHeight() + 5, 30, 30);
		this.add(coefficientList.get(0));
		compoundList.add(new Compound());
		Main.dataPanel.add(new DataPanel());
		Main.formulaPanel.add(new FormulaPanel());

		productBox.add(new Images("product.png", 846, 275));
		coefficientList.add(new JTextField(2));
		coefficientList.get(1).setBounds(881, 321, 30, 30);
		this.add(coefficientList.get(1));
		this.setLayout(null);
		compoundList.add(new Compound());
		Main.dataPanel.add(new DataPanel());
		Main.formulaPanel.add(new FormulaPanel());

		addReactant = new Images("addReactant.png", 32, 260);
		subReactant = new Images("subReactant.png", 32, 297);
		addProduct = new Images("addProduct.png", 1214, 260);
		subProduct = new Images("subProduct.png", 1214, 297);

		plus = new Images("plus.png", 0, 0);
	}

	// Resets all instance variables to their default values
	private void resetLists() {
		// TODO Auto-generated method stub
		compoundIndex = 0;
		reactantBox = new ArrayList<Images>();
		productBox = new ArrayList<Images>();
		coefficientList = new ArrayList<JTextField>();
		compoundList = new ArrayList<Compound>();
		fontFactor = 7;
	}

	// Displays the EquationPanel, which includes reactants, products, and a
	// console
	@Override
	public void paintComponent(Graphics g) {
		g.drawImage(background.img, background.x, background.y, null);
		int tempCompoundIndex = compoundIndex;
		for (int i = 0; i < reactantBox.size(); i++) {
			Font n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);
			compoundIndex = i;
			if (Main.formulaPanel.get(i).getLabel().length() > 0) {
				fontFactor = 585 / g.getFontMetrics().stringWidth(Main.formulaPanel.get(i).getLabel());
				fontFactor = fontFactor > 10 ? 10 : fontFactor;
			}
			n = new Font("SansSerif", Font.PLAIN, fontFactor * 2);
			g.setFont(n);
			g.drawImage(reactantBox.get(i).img, reactantBox.get(i).x, reactantBox.get(i).y, null);
			compoundList.get(i).paintComponent(g, reactantBox.get(i).x - reactantBox.get(i).img.getWidth() / 2 - 7,
					reactantBox.get(i).y + 25);
			if (reactantBox.size() > 1 && i < reactantBox.size() - 1) {
				g.drawImage(plus.img, reactantBox.get(i).x
						+ (reactantBox.get(i + 1).x + reactantBox.get(i + 1).img.getWidth() - reactantBox.get(i).x) / 2
						- plus.img.getWidth() / 2, reactantBox.get(i).y + 7, null);
			}
		}

		for (int i = 0; i < productBox.size(); i++) {
			Font n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);
			compoundIndex = i + reactantBox.size();
			if (Main.formulaPanel.get(compoundIndex).getLabel().length() > 0) {
				fontFactor = 585 / g.getFontMetrics().stringWidth(Main.formulaPanel.get(compoundIndex).getLabel());
				fontFactor = fontFactor > 10 ? 10 : fontFactor;
			}
			n = new Font("SansSerif", Font.PLAIN, fontFactor * 2);
			g.setFont(n);
			g.drawImage(productBox.get(i).img, productBox.get(i).x, productBox.get(i).y, null);
			compoundList.get(compoundIndex).paintComponent(g,
					productBox.get(i).x - productBox.get(i).img.getWidth() / 2 - 7, productBox.get(i).y + 25);
			if (productBox.size() > 1 && i < productBox.size() - 1) {
				g.drawImage(plus.img, productBox.get(i).x
						+ (productBox.get(i + 1).x + productBox.get(i + 1).img.getWidth() - productBox.get(i).x) / 2
						- plus.img.getWidth() / 2, productBox.get(i).y + 7, null);
			}
		}
		compoundIndex = tempCompoundIndex;
		g.drawImage(addReactant.img, addReactant.x, addReactant.y, null);
		g.drawImage(subReactant.img, subReactant.x, subReactant.y, null);
		g.drawImage(addProduct.img, addProduct.x, addProduct.y, null);
		g.drawImage(subProduct.img, subProduct.x, subProduct.y, null);
		Console.paintComponent(g);
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub

		int x = e.getX(), y = e.getY();

		// Add a reactant if the user presses the plus on the reactant side
		if (addReactant.pointOnImage(x, y) && reactantBox.size() < 4) {
			try {
				coefficientList.add(reactantBox.size(), new JTextField(2));
				reactantBox.add(new Images("reactant.png", 410 / (reactantBox.size() + 2), 275));
				for (int i = 0; i < reactantBox.size(); i++) {
					reactantBox.get(i).x = 640 / (reactantBox.size() + 1) * (i + 1) - 7;
					coefficientList.get(i).setBounds(reactantBox.get(i).x + reactantBox.get(i).img.getWidth() / 2 - 15,
							reactantBox.get(i).y + reactantBox.get(i).img.getHeight() + 5, 30, 30);
					this.add(coefficientList.get(i));
				}
				compoundList.add(reactantBox.size() - 1, new Compound());
				Main.dataPanel.add(reactantBox.size() - 1, new DataPanel());
				Main.formulaPanel.add(reactantBox.size() - 1, new FormulaPanel());
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} 
		
		// Subtract a reactant if the user presses the minus on the reactant side
		else if (subReactant.pointOnImage(x, y) && reactantBox.size() > 1) {
			this.remove(coefficientList.remove(reactantBox.size() - 1));
			compoundList.remove(reactantBox.size() - 1);
			Main.dataPanel.remove(reactantBox.size() - 1);
			Main.formulaPanel.remove(reactantBox.size() - 1);
			reactantBox.remove(reactantBox.size() - 1);
			for (int i = 0; i < reactantBox.size(); i++) {
				reactantBox.get(i).x = 640 / (reactantBox.size() + 1) * (i + 1) - 7;
				coefficientList.get(i).setBounds(reactantBox.get(i).x + reactantBox.get(i).img.getWidth() / 2 - 15,
						reactantBox.get(i).y + reactantBox.get(i).img.getHeight() + 5, 30, 30);
				this.add(coefficientList.get(i));
			}
		} 
		
		// Add a product if the user presses the plus on the product side
		else if (addProduct.pointOnImage(x, y) && productBox.size() < 4) {
			try {
				coefficientList.add(new JTextField(2));
				productBox.add(new Images("product.png", 0, 275));
				for (int i = productBox.size() - 1; i >= 0; i--) {
					productBox.get(i).x = 1280 - (640 / (productBox.size() + 1) * (i + 1) + 112);
					coefficientList.get(i + reactantBox.size()).setBounds(
							productBox.get(i).x + productBox.get(i).img.getWidth() / 2 - 15,
							productBox.get(i).y + productBox.get(i).img.getHeight() + 5, 30, 30);
					this.add(coefficientList.get(reactantBox.size() + productBox.size() - 1));

				}
				compoundList.add(new Compound());
				Main.dataPanel.add(new DataPanel());
				Main.formulaPanel.add(new FormulaPanel());
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		} 
		
		// Subtract a product if the user presses the minus on the product side
		else if (subProduct.pointOnImage(x, y) && productBox.size() > 1) {
			this.remove(coefficientList.get(productBox.size() - 1 + reactantBox.size()));
			coefficientList.remove(productBox.size() - 1 + reactantBox.size());
			compoundList.remove(productBox.size() - 1 + reactantBox.size());
			Main.dataPanel.remove(productBox.size() - 1 + reactantBox.size());
			Main.formulaPanel.remove(productBox.size() - 1 + reactantBox.size());
			productBox.remove(productBox.size() - 1);
			for (int i = productBox.size() - 1; i >= 0; i--) {
				productBox.get(i).x = 1280 - (640 / (productBox.size() + 1) * (i + 1) + 112);
				coefficientList.get(i + reactantBox.size()).setBounds(
						productBox.get(i).x + productBox.get(i).img.getWidth() / 2 - 15,
						productBox.get(i).y + productBox.get(i).img.getHeight() + 5, 30, 30);
			}

			// if the user clicked on the calculate button, calculate the
			// equation
		} 
		
		// If the user clicked on the calculate button, calculate the
		// equation
		else if (x > 534 && x < 721 && y > 388 && y < 435) {
			ArrayList<Compound> reactants = new ArrayList<Compound>();
			ArrayList<Compound> products = new ArrayList<Compound>();

			for (int i = 0; i < reactantBox.size(); i++)
				reactants.add(compoundList.get(i));

			for (int i = reactantBox.size(); i < compoundList.size(); i++)
				products.add(compoundList.get(i));

			ChemicalEquation equation = new ChemicalEquation(reactants, products);
			if (equation.solveEquation() && setCoefficients())
				Main.appState = 3;
		}

		// If the user presses on a reactant box, go to that compound's formulaPanel
		for (int i = 0; i < reactantBox.size(); i++) {
			if (reactantBox.get(i).pointOnImage(x, y)) {
				compoundIndex = i;
				fontFactor = 11;
				Main.appState = 1;
			}
		}

		// If the user presses on a product box, go to that compound's formulaPanel
		for (int i = 0; i < productBox.size(); i++) {
			if (productBox.get(i).pointOnImage(x, y)) {
				compoundIndex = i + reactantBox.size();
				fontFactor = 11;
				Main.appState = 1;
			}
		}
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	// When a key is typed, sets the coefficient 
	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

		for (int i = 0; i < compoundList.size(); i++) {
			compoundList.get(i).setCoefficient(new Integer(Integer.parseInt(coefficientList.get(i).getText())));
		}
	}

	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
	}

	// Sets the coefficients to the number in their corresponding text boxes
	public boolean setCoefficients() {
		for (int i = 0; i < compoundList.size(); i++) {
			if (coefficientList.get(i).getText().length() > 0) {
				try {
					compoundList.get(i).setCoefficient(Integer.parseInt(coefficientList.get(i).getText()));
				} catch (NumberFormatException n) {
					coefficientList.get(i).setText("");
					return false;
				}
			}
		}
		return true;
	}
}