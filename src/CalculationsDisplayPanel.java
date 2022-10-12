import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;

import javax.swing.JPanel;

public class CalculationsDisplayPanel extends JPanel implements MouseListener, KeyListener {

	Images calcDisplayBackground;

	// Call the constructor, which will enable keylistener, mouselistener, and
	// instantiate the background
	CalculationsDisplayPanel() throws IOException {
		this.addMouseListener(this);
		this.addKeyListener(this);
		this.requestFocus();
		calcDisplayBackground = new Images("CalcDisplayPanel.png", 0, 0);
	}

	// Displays the resulting data after the user presses "calculate"
	@Override
	public void paintComponent(Graphics g) {

		// Display the background
		g.drawImage(calcDisplayBackground.img, calcDisplayBackground.x, calcDisplayBackground.y, null);

		// Display the reactant labels and their corresponding data below them
		for (int i = 0; i < EquationPanel.reactantBox.size(); i++) {

			Font n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);
			EquationPanel.compoundIndex = i;

			if (Main.formulaPanel.get(i).getLabel().length() > 0) {
				EquationPanel.fontFactor = 585 / g.getFontMetrics().stringWidth(Main.formulaPanel.get(i).getLabel());
				EquationPanel.fontFactor = EquationPanel.fontFactor > 10 ? 10 : EquationPanel.fontFactor;
			}
			n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
			g.setFont(n);

			EquationPanel.compoundList.get(i)
					.paintComponent(g, EquationPanel.reactantBox.get(i).x
							- EquationPanel.reactantBox.get(i).img.getWidth() / 2 - 7,
					EquationPanel.reactantBox.get(i).y + 25);

			if (EquationPanel.reactantBox.size() > 1 && i < EquationPanel.reactantBox.size() - 1) {
				g.drawImage(EquationPanel.plus.img,
						EquationPanel.reactantBox.get(i).x
								+ (EquationPanel.reactantBox.get(i + 1).x
										+ EquationPanel.reactantBox.get(i + 1).img.getWidth()
										- EquationPanel.reactantBox.get(i).x) / 2
								- EquationPanel.plus.img.getWidth() / 2,
						EquationPanel.reactantBox.get(i).y + 7, null);
			}

			n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);
			g.drawString("(" + EquationPanel.compoundList.get(i).getCoefficient().toString() + ")",
					EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2 - 7,
					EquationPanel.reactantBox.get(i).y + EquationPanel.reactantBox.get(i).img.getHeight() + 5);

			String mass = "" + EquationPanel.compoundList.get(i).getMass();
			mass = mass.length() > 6 ? mass.substring(0, 6) : mass;
			mass += " g";

			String molarMass = "" + EquationPanel.compoundList.get(i).getMolarMass();
			molarMass = molarMass.length() > 6 ? molarMass.substring(0, 6) : molarMass;
			molarMass += " g/mol";

			String moles = "" + EquationPanel.compoundList.get(i).getMoles();
			moles = moles.length() > 6 ? moles.substring(0, 6) : moles;
			moles += " mol";

			if (i == 0) {
				g.drawString("Mass:",
						EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
								- g.getFontMetrics().stringWidth(mass) / 2 - 80,
						367);
				g.drawString("Molar Mass:",
						EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
								- g.getFontMetrics().stringWidth(molarMass) / 2 - 104,
						397);
				g.drawString("Moles:",
						EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
								- g.getFontMetrics().stringWidth(moles) / 2 - 77,
						427);
			}

			g.drawString(mass, EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
					- g.getFontMetrics().stringWidth(mass) / 2, 367);

			g.drawString(molarMass,
					EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
							- g.getFontMetrics().stringWidth(molarMass) / 2,
					397);

			g.drawString(moles, EquationPanel.reactantBox.get(i).x + EquationPanel.reactantBox.get(i).img.getWidth() / 2
					- g.getFontMetrics().stringWidth(moles) / 2, 427);
		}

		// Display the product labels, and their corresponding data below them

		for (int i = 0; i < EquationPanel.productBox.size(); i++) {
			Font n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);
			EquationPanel.compoundIndex = i + EquationPanel.reactantBox.size();
			if (Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel().length() > 0) {
				EquationPanel.fontFactor = 585
						/ g.getFontMetrics().stringWidth(Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel());
				EquationPanel.fontFactor = EquationPanel.fontFactor > 10 ? 10 : EquationPanel.fontFactor;
			}
			n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
			g.setFont(n);
			EquationPanel.compoundList.get(EquationPanel.compoundIndex).paintComponent(g,
					EquationPanel.productBox.get(i).x - EquationPanel.productBox.get(i).img.getWidth() / 2 - 7,
					EquationPanel.productBox.get(i).y + 25);
			if (EquationPanel.productBox.size() > 1 && i < EquationPanel.productBox.size() - 1) {
				g.drawImage(EquationPanel.plus.img,
						EquationPanel.productBox.get(i).x
								+ (EquationPanel.productBox.get(i + 1).x
										+ EquationPanel.productBox.get(i + 1).img.getWidth()
										- EquationPanel.productBox.get(i).x) / 2
								- EquationPanel.plus.img.getWidth() / 2,
						EquationPanel.productBox.get(i).y + 7, null);
			}
			n = new Font("SansSerif", Font.PLAIN, 13);
			g.setFont(n);

			// Draw coefficient under the compound
			g.drawString(
					"(" + EquationPanel.compoundList.get(i + EquationPanel.reactantBox.size()).getCoefficient()
							.toString() + ")",
					EquationPanel.productBox.get(i).x + EquationPanel.productBox.get(i).img.getWidth() / 2 - 7,
					EquationPanel.productBox.get(i).y + EquationPanel.productBox.get(i).img.getHeight() + 5);

			String mass = "" + EquationPanel.compoundList.get(i + EquationPanel.reactantBox.size()).getMass();
			mass = mass.length() > 6 ? mass.substring(0, 6) : mass;
			mass += " g";

			String molarMass = "" + EquationPanel.compoundList.get(i + EquationPanel.reactantBox.size()).getMolarMass();
			molarMass = molarMass.length() > 6 ? molarMass.substring(0, 6) : molarMass;
			molarMass += " g/mol";

			String moles = "" + EquationPanel.compoundList.get(i + EquationPanel.reactantBox.size()).getMoles();
			moles = moles.length() > 6 ? moles.substring(0, 6) : moles;
			moles += " g/mol";

			g.drawString(mass, EquationPanel.productBox.get(i).x + EquationPanel.productBox.get(i).img.getWidth() / 2
					- g.getFontMetrics().stringWidth(mass) / 2, 367);

			g.drawString(molarMass,
					EquationPanel.productBox.get(i).x + EquationPanel.productBox.get(i).img.getWidth() / 2
							- g.getFontMetrics().stringWidth(molarMass) / 2,
					397);

			g.drawString(moles, EquationPanel.productBox.get(i).x + EquationPanel.productBox.get(i).img.getWidth() / 2
					- g.getFontMetrics().stringWidth(moles) / 2, 427);
		}

	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
		// If the "NEW" button is pressed, create a new equation (reset
		// everything)
		int x = e.getX(), y = e.getY();
		if (x > 1161 && x < 1246 && y > 214 && y < 261) {
			try {
				Main.newEquation();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
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

}
