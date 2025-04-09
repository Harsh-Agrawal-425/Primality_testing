import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.swing.SwingUtilities;

public class PrimeGeneratorApplet extends Applet implements ActionListener {
    TextField bitField, attemptsField, timeField, avgTimeField, attemptsFieldAKS, timeFieldAKS, avgTimeFieldAKS;
    TextArea resultArea, rsaArea;
    Button generateButton, checkPrimeButton;
    String lastGeneratedPrime = "";

    public void init() {
        setLayout(new BorderLayout());

        Font appletFont = new Font("SansSerif", Font.BOLD, 20);
        Font resultFont = new Font("Monospaced", Font.PLAIN, 20);

        // Upper panel
        Panel upperPanel = new Panel(new BorderLayout());

        Panel topPanel = new Panel();
        Label bitLabel = new Label("Enter Bit Size:");
        bitLabel.setFont(appletFont);
        topPanel.add(bitLabel);

        bitField = new TextField(10);
        bitField.setFont(appletFont);
        topPanel.add(bitField);

        generateButton = new Button("Generate Prime");
        generateButton.setFont(appletFont);
        generateButton.addActionListener(this);
        topPanel.add(generateButton);

        upperPanel.add(topPanel, BorderLayout.NORTH);

        resultArea = new TextArea(10, 80);
        resultArea.setFont(resultFont);
        resultArea.setEditable(false);
        upperPanel.add(resultArea, BorderLayout.CENTER);

        Panel bottomPanel = new Panel(new GridLayout(2, 6));
        bottomPanel.add(createLabeledField("MR Attempts:", attemptsField = new TextField(10), appletFont));
        bottomPanel.add(createLabeledField("MR Time (ms):", timeField = new TextField(10), appletFont));
        bottomPanel.add(createLabeledField("MR Avg Time (ms):", avgTimeField = new TextField(10), appletFont));

        bottomPanel.add(createLabeledField("AKS Attempts:", attemptsFieldAKS = new TextField(10), appletFont));
        bottomPanel.add(createLabeledField("AKS Time (ms):", timeFieldAKS = new TextField(10), appletFont));
        bottomPanel.add(createLabeledField("AKS Avg Time (ms):", avgTimeFieldAKS = new TextField(10), appletFont));

        upperPanel.add(bottomPanel, BorderLayout.SOUTH);

        Panel checkPrimePanel = new Panel(new FlowLayout());
        checkPrimeButton = new Button("Check if Prime");
        checkPrimeButton.setFont(appletFont);
        checkPrimeButton.addActionListener(this);
        checkPrimePanel.add(checkPrimeButton);
        upperPanel.add(checkPrimePanel, BorderLayout.EAST);

        Panel lowerPanel = new Panel(new BorderLayout());
        rsaArea = new TextArea(10, 80);
        rsaArea.setFont(resultFont);
        rsaArea.setEditable(false);
        lowerPanel.add(new Label("RSA Key Generation performace for Miller Rabin and AKS:                                                                                                                                    Group 5(Arshit, Harsh and Prem)", Label.CENTER), BorderLayout.NORTH);
        lowerPanel.add(rsaArea, BorderLayout.CENTER);

        add(upperPanel, BorderLayout.NORTH);
        add(lowerPanel, BorderLayout.CENTER);
    }

    private Panel createLabeledField(String labelText, TextField textField, Font font) {
        Panel panel = new Panel();
        Label label = new Label(labelText);
        label.setFont(font);
        textField.setFont(font);
        textField.setEditable(false);
        panel.add(label);
        panel.add(textField);
        return panel;
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == generateButton) {
            resultArea.setText("");
            rsaArea.setText("");
            attemptsField.setText("");
            timeField.setText("");
            avgTimeField.setText("");
            attemptsFieldAKS.setText("");
            timeFieldAKS.setText("");
            avgTimeFieldAKS.setText("");

            generateButton.setEnabled(false);
            resultArea.setText("Generating prime... Please wait.\n\n");

            new Thread(new Runnable() {
                public void run() {
                    try {
                        int bits = Integer.parseInt(bitField.getText());
                        ProcessBuilder pb = new ProcessBuilder("./Prime_Generation_and_RSA/prime_generator", String.valueOf(bits));
                        pb.redirectErrorStream(true);
                        Process process = pb.start();

                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        final StringBuilder rsaOutput = new StringBuilder();
                        String line;
                        int lineCount = 0;

                        while ((line = reader.readLine()) != null) {
                            final String currentLine = line;
                            final int currentIndex = lineCount;

                            SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    switch (currentIndex) {
                                        case 0:
                                            resultArea.append("Miller-Rabin Prime: " + currentLine + "\n");
                                            lastGeneratedPrime = currentLine;
                                            break;
                                        case 1:
                                            attemptsField.setText(currentLine);
                                            break;
                                        case 2:
                                            timeField.setText(currentLine);
                                            break;
                                        case 3:
                                            avgTimeField.setText(currentLine);
                                            break;
                                        case 4:
                                            resultArea.append("\nAKS Prime: " + currentLine + "\n");
                                            break;
                                        case 5:
                                            attemptsFieldAKS.setText(currentLine);
                                            break;
                                        case 6:
                                            timeFieldAKS.setText(currentLine);
                                            break;
                                        case 7:
                                            avgTimeFieldAKS.setText(currentLine);
                                            break;
                                        default:
                                            rsaOutput.append(currentLine).append("\n");
                                            break;
                                    }
                                }
                            });

                            lineCount++;
                            try {
                                Thread.sleep(30);
                            } catch (InterruptedException ex) {
                                ex.printStackTrace();
                            }
                        }

                        SwingUtilities.invokeLater(new Runnable() {
                            public void run() {
                                rsaArea.setText(rsaOutput.toString());
                            }
                        });

                    } catch (final Exception ex) {
                        SwingUtilities.invokeLater(new Runnable() {
                            public void run() {
                                resultArea.setText("Execution error: " + ex.getMessage());
                            }
                        });
                    } finally {
                        SwingUtilities.invokeLater(new Runnable() {
                            public void run() {
                                generateButton.setEnabled(true);
                            }
                        });
                    }
                }
            }).start();
        }

        if (e.getSource() == checkPrimeButton) {
            checkPrimeButton.setEnabled(false);
            resultArea.append("\nChecking primality with AKS... Please wait.");

            new Thread(new Runnable() {
                public void run() {
                    try {
                        ProcessBuilder pb = new ProcessBuilder("./Prime_Validate/check_prime_with_aks");
                        pb.redirectErrorStream(true);
                        Process process = pb.start();

                        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
                        writer.write(lastGeneratedPrime + "\n");
                        writer.flush();
                        writer.close();

                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        final String result = reader.readLine();

                        SwingUtilities.invokeLater(new Runnable() {
                            public void run() {
                                resultArea.append("\nAKS Verification: " + lastGeneratedPrime + " is " +
                                        ("1".equals(result) ? "Prime" : "Composite"));
                                checkPrimeButton.setEnabled(true);
                            }
                        });

                    } catch (final Exception ex) {
                        SwingUtilities.invokeLater(new Runnable() {
                            public void run() {
                                resultArea.append("\nAKS Execution error: " + ex.getMessage());
                            }
                        });
                    }
                }
            }).start();
        }
    }
}
