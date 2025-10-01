import React from "react";
import { BrowserRouter as Router, Route, Routes, Link } from "react-router-dom";
import Home from "./pages/home";
import GenerateTrigger from "./pages/gen_trigger";
import CategorizeTrigger from "./pages/catag_trigger";
import { Container, Navbar, Nav } from "react-bootstrap";

const App = () => {
  return (
    <Router>
      {/* Bootstrap Navbar */}
      <Navbar bg="dark" variant="dark" expand="lg">
        <Container>
          <Navbar.Brand as={Link} to="/">DNA Uno</Navbar.Brand>
          <Navbar.Toggle aria-controls="basic-navbar-nav" />
          <Navbar.Collapse id="basic-navbar-nav">
            <Nav className="ms-auto">  {/* Aligns to the right */}
              <Nav.Link as={Link} to="/">Home</Nav.Link>
              <Nav.Link as={Link} to="/generate">Generate Template</Nav.Link>
              <Nav.Link as={Link} to="/categorize">Test Template</Nav.Link>
            </Nav>
          </Navbar.Collapse>
        </Container>
      </Navbar>

      {/* Main Content Area */}
      <Container className="mt-4">
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/generate" element={<GenerateTrigger />} />
          <Route path="/categorize" element={<CategorizeTrigger />} />
        </Routes>
      </Container>
    </Router>
  );
};

 export default App;
