import React, {useState} from "react";
import {Container, Form, Button, Table, Card} from "react-bootstrap";

const CategorizeTrigger = () => {
  const [inputSequences, setInputSequences] = useState("");
  const [sequenceType, setSequenceType] = useState("trigger"); // New state for type selection
  const [results, setResults] = useState([]);
  const [error, setError] = useState("");

  // Function to calculate GC content
  const calculateGCContent = (sequence) => {
    const gCount = (sequence.match(/G/g) || []).length;
    const cCount = (sequence.match(/C/g) || []).length;
    const totalLength = sequence.length;
    const gcContent = ((gCount + cCount) / totalLength) * 100;
    return gcContent.toFixed(2);
  };

  // Function to generate a random sequence of bases (A, T, C, G)
  const generateRandomBases = (length) => {
    const bases = ["A", "T", "C", "G"];
    return Array.from({length}, () => bases[Math.floor(Math.random() * bases.length)]).join("");
  };

  // Function to generate reverse complement of a DNA sequence
  const getReverseComplement = (sequence) => {
    // Create an object to map bases to their complements
    const complementMap = {
      A: "T",
      T: "A",
      C: "G",
      G: "C"
    };

    // Reverse the sequence and map each base to its complement
    const reverseComplement = sequence
        .split("")                      // Split sequence into an array of bases
        .reverse()                       // Reverse the array
        .map(base => complementMap[base]) // Map each base to its complement
        .join("");                       // Join the array back into a string
    return reverseComplement;
  };

  // // Function to generate complementary DNA sequence (A ↔ T, C ↔ G)
  // const getComplementarySequence = (sequence) => {
  //   return sequence
  //     .split("")
  //     .map((base) => {
  //       switch (base) {
  //         case "A": return "T";
  //         case "T": return "A";
  //         case "C": return "G";
  //         case "G": return "C";
  //         default: return base;
  //       }
  //     })
  //     .join("");
  // };

  const handleSubmit = (e) => {
    e.preventDefault();
    setError(""); // Clear errors

    // Split input into sequences
    const sequences = inputSequences
        .split(",")
        .map(seq => seq.trim())
        .filter(seq => seq.length > 0);

    if (sequences.length === 0) {
      setError("Please enter at least one valid DNA sequence.");
      return;
    }
    if (sequences.length > 5) {
      setError("You can enter a maximum of 5 sequences.");
      return;
    }

    // Validate sequences
    const generatedResults = sequences.map((seq, index) => {
      let templateSequence = seq;
      let isTrigger = sequenceType === "trigger";

      if (isTrigger) {
        // Ensure trigger is 10-15 bp
        if (seq.length < 10 || seq.length > 15 || !/^[CTAG]+$/.test(seq)) {
          setError(`Trigger sequence ${index + 1} must be 10-15 bases long and contain only C, T, A, and G.`);
          return null;
        }

        // Generate reverse complement for trigger sequence
        const reverseComplement = getReverseComplement(seq);

        // Display reverse complement in generated template
        templateSequence = reverseComplement; // Directly use reverse complement
        console.log(`Reverse complement generated: ${reverseComplement}`);
      } else {
        // Validate full template format
        const regex = /^([CTAG]{10,15})([CTAG]{4})(GACTC)(T)([CTAG]{10,15})$/;
        if (!regex.test(seq)) {
          setError(`Template sequence ${index + 1} is not in the correct format.`);
          return null;
        }
      }

      // Generate stats
      // const specificity = Math.floor(Math.random() * 30) + 70;
      const specificity = Math.floor(Math.random() * 30) + 70;
      // const Tm = (Math.random() * 10 + 60).toFixed(2);
      const Tm = (Math.random() * 10 + 50).toFixed(2);
      const gcContent = calculateGCContent(templateSequence);
      let status = gcContent < 40 || gcContent > 60 || Tm < 58 || Tm > 65 || specificity < 70 ? "Bad" : "Good";

      return {
        id: index + 1,
        originalSequence: seq,
        generatedTemplate: isTrigger ? templateSequence : "-", // Show generated template if trigger
        sequenceType: isTrigger ? "Trigger" : "Full Template",
        specificity: specificity + "%",
        Tm: Tm + "°C",
        gcContent: gcContent + "%",
        status,
      };
    }).filter(result => result !== null); // Remove null values from errors

    if (generatedResults.length > 0) {
      setResults(generatedResults.sort((a, b) => parseInt(b.specificity) - parseInt(a.specificity)));
    }
  };

  return (
      <Container className="mt-5 d-flex flex-column align-items-center">
        <h2 className="mb-4">Test Templates Sequences</h2>
        <p>Enter a <strong style={{color: "darkblue"}}>trigger (10-15 bp)</strong></p>
        <p style={{fontWeight: "bold", textAlign: "center"}}>or</p>
        <p>Enter a <strong style={{color: "darkblue"}}>full template</strong> directly in the following format:</p>

        <p>
          <strong>[</strong>10-15 trigger bases<strong>]</strong>
          <strong>[</strong>4 random bases<strong>]</strong>
          <strong>[</strong>GACTC (for cutting)<strong>]</strong>
          <strong>[</strong>T<strong>]</strong>
          <strong>[</strong>repeat trigger bases<strong>]</strong>
        </p>

        <br></br>
        <p>For testing: GTCGACTAATGCCAGACTCTGTCGACTAAT</p>

        {/* Error Message */}
        {error && <div className="alert alert-danger">{error}</div>}

        {/* Input Form */}
        <Form onSubmit={handleSubmit} className="w-50 p-4 border rounded shadow-sm bg-light">
          <Form.Group className="mb-3">
            <Form.Label>Choose Input Type:</Form.Label>
            <div className="d-flex">
              <Form.Check
                  type="radio"
                  label="Trigger (10-15 bp)"
                  name="sequenceType"
                  value="trigger"
                  checked={sequenceType === "trigger"}
                  onChange={() => setSequenceType("trigger")}
                  className="me-3"
              />
              <Form.Check
                  type="radio"
                  label="Full Template"
                  name="sequenceType"
                  value="template"
                  checked={sequenceType === "template"}
                  onChange={() => setSequenceType("template")}
              />
            </div>
          </Form.Group>

          <Form.Group className="mb-3">
            <Form.Label>Enter Up to 5 Sequences (comma-separated):</Form.Label>
            <Form.Control
                type="text"
                value={inputSequences}
                onChange={(e) => setInputSequences(e.target.value)}
                required
                placeholder="Example: ATGCGATCGATC"
            />
          </Form.Group>
          <Button variant="primary" type="submit" className="w-100">Process</Button>
        </Form>

        {/* Results Table */}
        {results.length > 0 && (
            <Card className="mt-4 w-75 shadow-sm" style={{maxHeight: "60vh", overflowY: "auto"}}>
              <Card.Body>
                <Card.Title>Sequence Analysis Results</Card.Title>
                <div style={{overflowX: "auto"}}>
                  <Table striped bordered hover className="table-responsive">
                    <thead>
                    <tr>
                      <th>Rank</th>
                      <th>Input Type</th>
                      <th style={{minWidth: "200px"}}>Original Sequence</th>
                      <th style={{minWidth: "250px"}}>Generated Template</th>
                      <th>Specificity</th>
                      <th>Tm (°C)</th>
                      <th>GC Content</th>
                      <th>Status</th>
                    </tr>
                    </thead>
                    <tbody>
                    {results.map((result, index) => (
                        <tr key={result.id}>
                          <td>{index + 1}</td>
                          <td>{result.sequenceType}</td>
                          <td style={{wordBreak: "break-word", whiteSpace: "normal"}}>{result.originalSequence}</td>
                          <td style={{wordBreak: "break-word", whiteSpace: "normal"}}>Trigger reverse
                            complement: {result.generatedTemplate}</td>
                          <td>{result.specificity}</td>
                          <td>{result.Tm}</td>
                          <td>{result.gcContent}</td>
                          <td>{result.status}</td>
                        </tr>
                    ))}
                    </tbody>
                  </Table>
                </div>
              </Card.Body>
            </Card>
        )}
        {/* White Section at the Bottom */}
        <div style={{ height: "100px", backgroundColor: "white", width: "100%" }}></div>
      </Container>
  );
};


export default CategorizeTrigger;
