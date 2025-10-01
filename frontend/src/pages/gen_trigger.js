import React, { useState } from "react";
import {Container, Form, Button, Card, Tooltip, OverlayTrigger} from "react-bootstrap";

const GenerateTrigger = () => {
  const [dnaSequence, setDnaSequence] = useState("");
  const [temperature, setTemperature] = useState("");
  const [generatedTrigger, setGeneratedTrigger] = useState("");
  const [fullSequence, setFullSequence] = useState("");
  const [tm1, setTm1] = useState("");
  const [numBonds, setNumBonds] = useState("");
  const [error, setError] = useState(""); // Added state for error message

  // Function to generate a random sequence of bases (A, T, C, G)
  const generateRandomBases = (length) => {
    const bases = ["A", "T", "C", "G"];
    return Array.from({ length }, () => bases[Math.floor(Math.random() * bases.length)]).join("");
  };

  const handleSubmit = (e) => {
    e.preventDefault();

        // Check each sequence for valid characters (CTAG)
    if (!/^[CTAG]+$/.test(dnaSequence)) {
      setError("DNA sequence contains invalid characters. Only 'C', 'T', 'A', and 'G' are allowed.");
      return;
    }

    if (dnaSequence.length < 20) {
      alert("DNA sequence must be at least 20 bases long.");
      return;
    }

    const randomBases = generateRandomBases(4); // Generate 4 random bases

    // Simulating generated results
    setGeneratedTrigger(dnaSequence.slice(10, 22)); // Example logic, just cutting from pos 11- 22
    setFullSequence(`${dnaSequence.slice(0, 10)}-${randomBases}-GACTC-T-${dnaSequence.slice(0, 10)}`);
    // setTm1(Math.random() * 30 + 50); // Example random TM1
    setTm1(Math.random() * 10 + 50);
    setNumBonds(Math.floor(Math.random() * 10 + 2)); // Example number of bonds
    setError(""); // Clear error if input is valid
  };

  return (
    <Container className="mt-5 d-flex flex-column align-items-center">
      <h2 className="mb-4">Generate Template</h2>
      <p> Please input your DNA sequence below. Minimum sequence length must be greater than 20 bases.</p>
      <Form onSubmit={handleSubmit} className="w-50 p-4 border rounded shadow-sm bg-light">
        <Form.Group className="mb-3">
          <Form.Label>Enter DNA Sequence:</Form.Label>
          <Form.Control
            type="text"
            value={dnaSequence}
            onChange={(e) => setDnaSequence(e.target.value)}
            required
            placeholder="Enter DNA sequence here... 5'-3'"
          />
          {error && <p className="text-danger">{error}</p>} {/* Display error message */}
        </Form.Group>

        <Form.Group className="mb-3">
          <Form.Label>Temperature of Reaction (°C):</Form.Label>
          <Form.Control
            type="number"
            value={temperature}
            onChange={(e) => setTemperature(e.target.value)}
            required
            placeholder="Enter temperature in °C"
          />
        </Form.Group>

        <Button variant="success" type="submit" className="w-100">
          Generate
        </Button>
      </Form>

      {/* Display Generated Results */}
      {generatedTrigger && (
        <Card className="mt-4 w-50 shadow-sm">
          <Card.Body>
            <Card.Title>Results</Card.Title>
            <Card.Text><strong>Generated Trigger (10-15bp):</strong> {generatedTrigger}</Card.Text>
            <Card.Text><strong>Full Sequence:</strong> {fullSequence}
            <OverlayTrigger
                placement="right-end"
                overlay={
                  <Tooltip>Understanding the generated format: [10-15 trigger bases] - [4 random bases] -
                    [GACTC (for cutting)] - [T] - [repeat trigger bases]
                  </Tooltip>
                }
            >
              <span
                  style={{
                    display: "inline-flex",
                    alignItems: "center",
                    justifyContent: "center",
                    width: "18px",
                    height: "18px",
                    borderRadius: "50%",
                    backgroundColor: "lightgrey",
                    color: "black",
                    fontSize: "12px",
                    fontWeight: "bold",
                    cursor: "pointer",
                    marginLeft: "5px"
                  }}
              >
            i
          </span>
            </OverlayTrigger>
            </Card.Text>

            <Card.Text><strong>TM1 (°C):</strong> {tm1.toFixed(2)}</Card.Text>
            <Card.Text><strong># Bonds:</strong> {numBonds}</Card.Text>
          </Card.Body>
        </Card>
      )}
    </Container>
  );
};

export default GenerateTrigger;
