import java.io.IOException;
import java.util.ArrayList;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class Main {
	static CalculationsDisplayPanel calcDisplayPanel; 
	static EquationPanel equationPanel;
	static ArrayList<DataPanel> dataPanel;
	static ArrayList<FormulaPanel> formulaPanel;
	static int appState = 0;
	// NOTE: May need to change the existence of "NEW" for the gamestates,
	// particularly for equation, also maybe make a method to initialize the
	// panels instead of doing it in each switch

	public static void main(String[] args) throws IOException, InterruptedException {
		// Declare and initialize a JFrame
		JFrame myFrame = new JFrame();
		
		// Make the frame visible
		myFrame.setTitle("StoichBuddy");
		myFrame.setSize(1280, 730);
		myFrame.setResizable(false);
		myFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		myFrame.setLocationRelativeTo(null);
		myFrame.setVisible(true);

		// Declare all panels
		newEquation();
		
		// Create a panel
		JPanel myPanel = new JPanel();
		myPanel = equationPanel;

		myPanel.setSize(1280, 740);
		myFrame.add(myPanel);
		myPanel.setVisible(true);
		
		while (true) {
			
			switch (appState) {
			case 0:
				// Set panel to Equation Panel
				myPanel = equationPanel;
				break;
			case 1:
				// Set panel to FormulaPnel
				myPanel = formulaPanel.get(EquationPanel.compoundIndex);

				break;
			case 2:
				// Set panel to Data Panel
				myPanel = dataPanel.get(EquationPanel.compoundIndex);
				break;
			case 3: 
				// Set panel to CalculationsDisplayPanel
				myPanel = calcDisplayPanel;
				break;
			}

			// Put the panel inside the frame
			myFrame.add(myPanel);
			myFrame.revalidate();

			
			myPanel.setVisible(true);
			myPanel.requestFocus();

			// Process in state, until switch
			int cur = appState;

			while (appState == cur) {
				Thread.sleep(60);
				myPanel.repaint();
			}
			
			myPanel.setVisible(false);
			myFrame.remove(myPanel);
		}
	}

	public static void newEquation() throws IOException {
		// TODO Auto-generated method stub
		dataPanel = new ArrayList<DataPanel>();
		DataPanel.dataBackground = new Images("DataPanel.png", 0, 0);
		formulaPanel = new ArrayList<FormulaPanel>();
		FormulaPanel.formulaBackground = new Images("FormulaPanel.png", 0, 0);
		equationPanel = new EquationPanel();
		calcDisplayPanel = new CalculationsDisplayPanel(); 
		appState = 0;
	}
}