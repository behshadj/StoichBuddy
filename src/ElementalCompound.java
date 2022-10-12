import java.awt.Font;
import java.awt.Graphics;
import java.io.IOException;

public class ElementalCompound extends Element {
	private boolean completed; // Checks if end brace has been entered for input
	private double molarMass;
	private Images blankSubscript; 
	
	
	public ElementalCompound() throws IOException {
		super(0, 0);
		blankSubscript = new Images("subscriptBox.png", 0, 0);
	}

	// Paints elemental compound based on input
	@Override
	public void paintComponent(Graphics g, int i, int width, int y) {
		width += g.getFontMetrics().stringWidth("(");

		for (int j = 0; j < getElement().size(); j++) {
			Font n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
			g.setFont(n);
			width += g.getFontMetrics().stringWidth(getElement().get(j).getAbbr() + " ");
			 n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor);
			g.setFont(n);
			getElement().get(j).paintComponent(g, i, width, y);

		}
		Font n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
		g.setFont(n);
		width += g.getFontMetrics().stringWidth(") ");
		 n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor);
		g.setFont(n);
		if (this.getSubscript() > 1)
			g.drawString(this.getSubscript().toString(), width - 6, y + 7);
		else if (this.getSubscript() == 0 && isCompleted()) {
			g.drawImage(this.blankSubscript.img, width - 6, y, null);
		}
		n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
		g.setFont(n); 
	}

	//get molar mass of elemental compound by summing the elements
	@Override
	public double getMolarMass() {
		this.molarMass = 0;
		for (int i = 0; i < getElement().size(); i++) {
			this.molarMass += getElement().get(i).getMolarMass();
		}
		return this.molarMass * this.getSubscript();
	}

	//getters and setters
	public boolean isCompleted() {
		return completed;
	}

	public void setCompleted(boolean completed) {
		this.completed = completed;
	}
}

