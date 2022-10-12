import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;

public class Compound {

	private Integer coefficient;
	private ArrayList<Element> elementList = new ArrayList<Element>();
	private double molarMass, moles, density, volume, mass;

	// Call the constructor, which sets the coefficient to the default value of
	// 1
	Compound() {
		this.setCoefficient(new Integer(1));
	}

	// Display the compound (i.e. Element abbreviations, subscripts)
	public void paintComponent(Graphics g, int x, int y) {
		int labelX = (210
				- g.getFontMetrics().stringWidth(Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel())) / 2
				+ x;
		Font n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
		g.setFont(n);
		g.drawString(Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel(), labelX, y);
		
		for (int i = 0; i < EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.size(); i++) {
			String add = "";
			if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.get(i).isSingle())
				add = Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel().substring(0,
						EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.get(i).getIndex()
								+ EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.get(i)
										.getAbbr().length()
								+ 1);
			else {
				add = Main.formulaPanel.get(EquationPanel.compoundIndex).getLabel().substring(0,
						EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.get(i).getIndex());
			}

			EquationPanel.compoundList.get(EquationPanel.compoundIndex).elementList.get(i).paintComponent(g, i,
					labelX + g.getFontMetrics().stringWidth(add), y);
		}
	}

	// Recalculate the molar mass using the current set of elements
	public void setMolarMass() {
		this.molarMass = 0;
		for (int i = 0; i < elementList.size(); i++) {
			this.molarMass += elementList.get(i).getMolarMass();
		}
	}

	// Accessors and mutators of all instance variables in this class 
	public Integer getCoefficient() {
		return coefficient;
	}

	public void setCoefficient(Integer coefficient) {
		this.coefficient = coefficient;
	}

	public ArrayList<Element> getElementList() {
		return elementList;
	}

	public double getMolarMass() {
		return molarMass;
	}

	public double getVolume() {
		return volume;
	}

	public void setVolume(double volume) {
		this.volume = volume;
	}

	public double getMoles() {
		return moles;
	}

	public void setMoles(double moles) {
		this.moles = moles;
	}

	public double getMass() {
		return mass;
	}

	public void setMass(double mass) {
		this.mass = mass;
	}

	public double getDensity() {
		return density;
	}

	public void setDensity(double density) {
		this.density = density;
	}
}
