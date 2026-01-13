Programs for aerospace applications.

Split into aerodynamics, stability, and propulsion categories (control code may be added, but probably not cause it's harder than using simulink).

- Aerodynamics contains incompressible (potential flow theory) and compressible applications.
  - also contains classes for curve design. useful for wing sections and creating geometries for CFD domains
- Propulsion contains programs that automate the jet engine design process.
  - the engine class is intended to wrap the preliminary design steps into a simply-to-manipulate object
    - flight performance (constraint & mission analysis)
    - design point analysis
    - off-design analysis
- Stability contains a calculator for stability derivatives based on aircraft geometry
  - the geometry must be provided first and each value input into the program
  - the process is iterative, but the stability coefficients can be recalculated quickly; the geometry being the only thing laborious to iterate (potentially made easier with Onshape variable studios)
  - lacks a complementary control program
