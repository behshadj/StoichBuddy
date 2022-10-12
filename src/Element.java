import java.awt.Graphics;
import java.util.ArrayList;

public abstract class Element {
	private int id;
	private Integer subscript = new Integer(1);
	private String abbr;
	private double molarMass;
	private int index;
	private boolean single;
	private ArrayList<Element> element = new ArrayList<Element>();

	// Calls the constructor, which instantiates the element's Atomic Number,
	// subscript, Abbreviation, and Molar Mass
	Element(int i, int s) {
		this.setId(i);
		this.setSubscript(s);
		this.setAbbr(AtomicData.abbreviation[i]);
		this.setMolarMass(AtomicData.mMass[i] * s);
	}

	// Once implemented, will display the Element's label and substring
	public abstract void paintComponent(Graphics g, int i, int width, int y);

	// All accessors and mutators for Element
	public abstract double getMolarMass();

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public Integer getSubscript() {
		return subscript;
	}

	public void setSubscript(Integer subscript) {
		this.subscript = subscript;
	}

	public String getAbbr() {
		return abbr;
	}

	public void setAbbr(String abbr) {
		this.abbr = abbr;
	}

	public void setMolarMass(double molarMass) {
		this.molarMass = molarMass;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public boolean isSingle() {
		return single;
	}

	public void setSingle(boolean single) {
		this.single = single;
	}

	public ArrayList<Element> getElement() {
		return element;
	}

	public void setElement(ArrayList<Element> element) {
		this.element = element;
	}

}