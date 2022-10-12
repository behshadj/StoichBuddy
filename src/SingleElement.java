import java.awt.Font;
import java.awt.Graphics;
import java.io.IOException;

public class SingleElement extends Element {

	Images blankSubscript;
	
	// Calls the constructor, which instantiates the instance variables 
	SingleElement(int i, int s) throws IOException {
		super(i, s);
		blankSubscript = new Images("subscriptBox.png", 0, 0);
	}

	// Returns the abbreviation
	public String toString() {
		return getAbbr();
	}

	// Displays the Element and its subscript
	@Override
	public void paintComponent(Graphics g, int i, int width, int y) {

		if (getSubscript() == 0) {
			g.drawImage(blankSubscript.img, width - 4, y, null);
		} else if (getSubscript() > 1) {
			Font n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor);
			g.setFont(n);
			g.drawString(getSubscript().toString(), width - 6, 7 + y);
			n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
			g.setFont(n);
		}
	}

	// Return the molar mass
	@Override
	public double getMolarMass() {
		// TODO Auto-generated method stub
		return AtomicData.mMass[this.getId()] * this.getSubscript();
	}
}
