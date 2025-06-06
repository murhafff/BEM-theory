import java.util.Arrays;

public class Main {
    // Linear interpolation helper function: given arrays x[] and y[], and a query xi,
    // returns linearly interpolated value of y at x = xi. If xi is outside the range of x[],
    // it clamps to the nearest endpoint of y[].
    static double interpolate(double[] x, double[] y, double xi) {
        for (int i = 0; i < x.length - 1; i++) {
            if (xi >= x[i] && xi <= x[i + 1]) {
                double t = (xi - x[i]) / (x[i + 1] - x[i]);
                return y[i] + t * (y[i + 1] - y[i]);
            }
        }
        if (xi < x[0]) {    // If xi is below x[0], return y[0]; if above x[last], return y[last].
            return y[0];
        } else {
            return y[y.length - 1];
        }
    }

    public static void main(String[] args) {

        // *******************************
        // ***** GIVEN CONSTANTS *********
        // *******************************
        final double rho = 0.4135;          // Air density rho at 10,000 m altitude (kg/m^3).
        final double u0 = 30.0;             // Freestream (flight) velocity u0 in m/s, as specified (30 m/s).
        final double aSound = 299.5;        // Speed of sound at 10,000 m altitude (approximate, in m/s).Used to check Mach number at the tip.
        final int segments = 100;           // Number of radial segments for BEM discretization
        final double rHubFraction = 0.1;    // Hub location as a fraction of radius: r_hub_fraction = 0.1 (10% of total radius).
        final double maxDiameter = 3.0;     // Maximum allowable diameter in meters.
        // Thrust target and allowable tolerance (100 N ± 1%).
        final double thrustTarget = 100.0;
        final double thrustLowerBound = thrustTarget * 0.99; // 99 N
        final double thrustUpperBound = thrustTarget * 1.01; // 101 N

        // *******************************
        // ***** AIRFOIL POLAR DATA ******
        // *******************************

        double[] AoA_table = {-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18};                                                       // Angles of attack (degrees) for the NACA 0012 polar.
        double[] CL_table  = {-0.2187, 0.0, 0.2188, 0.4366, 0.6528, 0.8596, 1.0715, 1.2671, 1.4461, 1.5998, 1.6895};        // Corresponding lift coefficients (CL) for NACA 0012 at those AoA values.
        double[] CD_table  = { 0.0096, 0.0102, 0.0114, 0.0132, 0.0158, 0.0193, 0.0241, 0.0307, 0.0390, 0.0553, 0.0960};     // Corresponding drag coefficients (CD) for NACA 0012 at those AoA values.

        // *********************************************************
        // ***** BRUTE-FORCE SEARCH OVER DESIGN VARIABLES **********
        // *********************************************************
        // Loop over discrete sets of:
        //   1) Number of blades B (integer values 2, 3, or 4).
        //   2) Propeller radius (meters) from 0.5 to 1.5 in steps of 0.1 (so diameter ≤ 3.0 m).
        //   3) Chord-scale factors from 0.5 to 1.5 in steps of 0.1, applied uniformly to a base chord distribution.
        //   4) Shaft speed RPM from 100 to 10000 in steps of 100 rpm.
        // For each combination:
        //   - Build the blade geometry (r[], chord[]).
        //   - Compute an optimized twist distribution θ(r) via the new twist-refinement procedure.
        //   - Check diameter constraint (2 * radius ≤ maxDiameter).
        //   - Compute tip Mach number and skip if Mach ≥ 0.8.
        //   - Run a final BEM analysis using the refined θ(r).
        //   - Compute total thrust T.
        //   - If T ∈ [95,105] N, we stop and report that design.

        boolean foundValidDesign = false;   // Flag indicating whether we have found a valid design.

        // Variables to store the successful design parameters for printing later.
        int bestB = -1;
        double bestRadius = -1; double bestChordScale = -1;
        int bestRPM = -1;
        double bestThrust = 0; double bestTorque = 0; double bestPower = 0; double bestEfficiency = 0; double J = 0;
        double[] chord = new double[segments];      // chord[i]
        double[] aValues = new double[segments]; // axial induction
        double[] a_omegaValues = new double[segments]; // tangential induction
        double[] PhiValues = new double[segments]; // advance angle
        double[] betaValues = new double[segments]; // blade-pitch distribution

        for (int B = 2; B <= 4 && !foundValidDesign; B++) {                                 // Begin looping over number of blades B: try 2, then 3, then 4.
            for (double radius = 0.6; radius <= 1.5 && !foundValidDesign; radius += 0.05) {  // Loop over radius from 0.5 to 1.5 meters in 0.1 m increments.
                double diameter = 2.0 * radius;                                             // Check diameter constraint: diameter = 2 * radius must not exceed maxDiameter.
                if (diameter > maxDiameter) {                                               // If diameter > 3.0 m, skip this radius.
                    continue;
                }

                double r_hub = rHubFraction * radius;       // Hub radius is 10% of the propeller radius.
                double dr = (radius - r_hub) / segments;    // the blade is discretized from r_hub to radius.
                double[] r = new double[segments];          // Build the radial coordinate array r[i] using the mid-point rule.
                for (int i = 0; i < segments; i++) {
                    r[i] = r_hub + dr * (i + 0.5);
                }

                for (double chordScale = 0.5; chordScale <= 2.5 && !foundValidDesign; chordScale += 0.05) {  // determine the size of the chord.
                    double[] chordBase = new double[segments];                                              // Pre-compute the base chord distribution for all stations.
                    for (int i = 0; i < segments; i++) {                                                    // Base chord distribution, tapering linearly from 0.4 m at inboard to 0.3 m at outboard.
                        chordBase[i] = 0.15 - (0.10 * i / segments);
                    }
                    for (int i = 0; i < segments; i++) {
                        chord[i] = chordBase[i] * chordScale;   // Then the actual chord[i] = chord_base[i] * chordScale.
                        if (chord[i] < 0.01) {                  // Ensure chord is never negative or zero; clamp to a small positive if needed.
                            chord[i] = 0.01;
                        }
                    }

                    for (int RPM = 100; RPM <= 10000 && !foundValidDesign; RPM += 25) {    // Loop over RPM from 100 to 10000 in steps of 100 rpm.
                        double omega = RPM * Math.PI / 30.0;                                // Convert RPM to angular speed omega in rad/s.

                        double tipSpeed = omega * radius;                                   // Compute Mach number at the tip.
                        double machTip = tipSpeed / aSound;
                        if (machTip > 0.8) {                                                // Skip if tip Mach number is ≥ 0.8.
                            continue;
                        }

                        // **************************
                        // **** TWIST-REFINEMENT ****
                        // **************************

                        double alphaOptDeg = AoA_table[0];
                        double maxLD = CL_table[0] / CD_table[0];
                        for (int i = 1; i < AoA_table.length; i++) {        // 1) Find α_optimum (radians): angle of attack that maximizes CL/CD.
                            double LD = CL_table[i] / CD_table[i];
                            if (LD > maxLD) {
                                maxLD = LD;
                                alphaOptDeg = AoA_table[i];
                            }
                        }
                        double alphaOpt = Math.toRadians(alphaOptDeg);
                        // 2) Allocate arrays for induction factors (a, a_omega), the inflow angle (Phi) and pitch angle (beta).

                        for (int i = 0; i < segments; i++) {                // 3) Initial guess: beta_i = Phi0_i + α_opt,  where Phi0_i = arctan(u0/(ω*r[i]))
                            double Phi0 = Math.atan(u0 / (omega * r[i]));
                            betaValues[i] = Phi0 + alphaOpt;
                        }

                        // 4) Twist-refinement loop: update θ until max|Δθ| < tolerance or maxRefines reached.
                        final double tolerancebeta = Math.toRadians(0.1);  // 0.1° tolerance
                        for (int refine = 0; refine < 200; refine++) {
                            for (int i = 0; i < segments; i++) {    // 4a) Run a BEM pass with the current beta(r):
                                double a  = 0.1;
                                double a_omega = 0.01;

                                for (int iter = 0; iter < 10000; iter++) {
                                    double Phi = Math.atan(((1 + a) * u0) / ((1 - a_omega) * omega * r[i]));
                                    double alphaLocal = betaValues[i] - Phi;
                                    double alphaDeg = Math.toDegrees(alphaLocal);

                                    double CL = interpolate(AoA_table, CL_table, alphaDeg);
                                    double CD = interpolate(AoA_table, CD_table, alphaDeg);

                                    double sinPhi = Math.sin(Phi);
                                    double cosPhi = Math.cos(Phi);

                                    double X = (B * chord[i] * (CL * cosPhi - CD * sinPhi)) / (8.0 * Math.PI * r[i] * sinPhi * sinPhi);// equation for a
                                    double Y = (B * chord[i] * (CL * sinPhi + CD * cosPhi)) / (8.0 * Math.PI * r[i] * sinPhi * cosPhi);// equation for a_omega

                                    double a_initial  = X / (1 - X);
                                    double aNew  = 0.5 * (a + a_initial);
                                    double a_omega_Initial  = Y / (1 + Y);
                                    double a_omegaNew = 0.5 * (a_omega + a_omega_Initial);


                                    if (Math.abs(aNew - a) < 1e-6 && Math.abs(a_omegaNew - a_omega) < 1e-6) {
                                        a  = aNew;
                                        a_omega = a_omegaNew;
                                        break;
                                    }
                                    a  = aNew;
                                    a_omega = a_omegaNew;
                                }

                                aValues[i]   = a;
                                a_omegaValues[i]  = a_omega;
                                PhiValues[i] = Math.atan2((1 + a) * u0, (1 - a_omega) * omega * r[i]);

                            }

                            double maxdifference = 0.0;                             // 4b) Update beta_new(i) = Phi(i) + α_opt and check convergence.
                            double[] newbeta = new double[segments];
                            for (int i = 0; i < segments; i++) {
                                newbeta[i] = betaValues[i]+ 0.5* (PhiValues[i] + alphaOpt - betaValues[i]);
                                double difference = Math.abs(newbeta[i] - betaValues[i]);
                                if (difference > maxdifference) {
                                    maxdifference = difference;
                                }

                            }

                            System.arraycopy(newbeta, 0, betaValues, 0, segments);  // Copy newbeta into beta and break if converged.
                            if (maxdifference < tolerancebeta) {
                                break;
                            }
                        }
                        boolean betaCheck = true;
                        double betaLimit = Math.toRadians(70);
                        for (int i = 0; i < segments; i++) {
                            if (betaValues[i] >= betaLimit) {
                                betaCheck = false;
                                break;
                            }
                        }
                        if (!betaCheck) {
                            System.out.println("over");
                            // If any station's final theta exceeds 70°, reject this combination
                            continue;  // move to the next RPM
                        }

                        // ***********************
                        // **** FINAL BEM RUN ****
                        // ***********************
                        // (Recompute aVals[], apVals[], phiVals[] one last time with converged beta(r))
                        for (int i = 0; i < segments; i++) {
                            double a  = 0.1;
                            double a_omega = 0.01;
                            for (int iter = 0; iter < 10000; iter++) {
                                double Phi = Math.atan(((1 + a) * u0) / ((1 - a_omega) * omega * r[i]));
                                double alphaLocal = betaValues[i] - Phi;
                                double alphaDeg = Math.toDegrees(alphaLocal);

                                double CL = interpolate(AoA_table, CL_table, alphaDeg);
                                double CD = interpolate(AoA_table, CD_table, alphaDeg);

                                double sinPhi = Math.sin(Phi);
                                double cosPhi = Math.cos(Phi);

                                double X = (B * chord[i] * (CL * cosPhi - CD * sinPhi)) / (8.0 * Math.PI * r[i] * (sinPhi*sinPhi));// equation for a
                                double Y = (B * chord[i] * (CL * sinPhi + CD * cosPhi)) / (8.0 * Math.PI * r[i] * sinPhi * cosPhi);// equation for a_omega

                                double a_initial  = X / (1 - X);
                                double aNew  = 0.5 * (a + a_initial);
                                double a_omega_Initial  = Y / (1 + Y);
                                double a_omegaNew = 0.5 * (a_omega + a_omega_Initial);

                                if (Math.abs(aNew - a) < 1e-6 && Math.abs(a_omegaNew - a_omega) < 1e-6) {
                                    a  = aNew;
                                    a_omega = a_omegaNew;
                                    break;
                                }
                                a  = aNew;
                                a_omega = a_omegaNew;
                            }
                            aValues[i]   = a;
                            a_omegaValues[i]  = a_omega;
                            PhiValues[i] = Math.atan(((1 + a) * u0) / ((1 - a_omega) * omega * r[i]));
                        }

                        // *******************************
                        // ***** THRUST & TORQUE *********
                        // *******************************
                        double Tsum = 0.0;
                        double Qsum = 0.0;
                        for (int i = 0; i < segments; i++) {
                            double alphaLocal = betaValues[i] - PhiValues[i];            // local AoA
                            double alphaDeg = Math.toDegrees(alphaLocal);      // convert to degrees
                            double CL = interpolate(AoA_table, CL_table, alphaDeg);
                            double CD = interpolate(AoA_table, CD_table, alphaDeg);

                            double sinPhi = Math.sin(PhiValues[i]);
                            double cosPhi = Math.cos(PhiValues[i]);

                            // "common" dynamic pressure * blade area factor:
                            // 0.5 * rho * u0^2 * B * chord[i] * (1 + a)^2 / sin^2(phi)
                            double common = 0.5 * rho * (u0*u0) * B * chord[i] * ((1 + aValues[i])*(1 + aValues[i]))
                                    / (sinPhi * sinPhi);

                            // Sectional thrust dT = common * (CL*cos(phi) - CD*sin(phi)) * dr
                            double dT = common * (CL * cosPhi - CD * sinPhi) * dr;
                            // Sectional torque dQ = common * (CL*sin(phi) + CD*cos(phi)) * r[i] * dr
                            double dQ = common * (CL * sinPhi + CD * cosPhi) * r[i] * dr;

                            // Accumulate sums
                            Tsum += dT;
                            Qsum += dQ;
                        }

                        // Compute required shaft power = Qsum * omega (in Watts).
                        double power = Qsum * omega;

                        // Propulsive efficiency η = (T * u0) / P.
                        double efficiency = (Tsum * u0) / power;

                        // Check if total thrust Tsum is within ±10% of 100 N.
                        if (Tsum >= thrustLowerBound && Tsum <= thrustUpperBound) {
                            // We found a valid design. Record its parameters.
                            bestB = B;
                            bestRadius = radius;
                            bestChordScale = chordScale;
                            bestRPM = RPM;
                            bestThrust = Tsum;
                            bestTorque = Qsum;
                            bestPower = power;
                            bestEfficiency = efficiency;
                            foundValidDesign = true;
                            J = u0/((RPM/60) * radius*2);
                            // Break out of the RPM loop.
                            break;
                        }
                    } // End RPM loop

                    // *** END OF REMOVED: constant-pitch loop ***
                    // } // End pitch loop (was here)
                } // End chord-scale loop
            } // End radius loop
        } // End B loop

        // *******************************
        // ***** OUTPUT RESULTS **********
        // *******************************

        if (foundValidDesign) {
            // Print the successful design parameters and performance.
            System.out.println("=== VALID PROPELLER DESIGN FOUND ===");
            System.out.println("a(r): " + Arrays.toString(aValues));
            System.out.println("a'(r): " + Arrays.toString(a_omegaValues));
            System.out.println("chord: " + Arrays.toString(chord));
            System.out.println("phi(r) [deg]: " + Arrays.toString(
                    Arrays.stream(PhiValues).map(p->Math.toDegrees(p)).toArray()));
            System.out.println("beta(r) [deg]: " + Arrays.toString(
                    Arrays.stream(betaValues).map(p->Math.toDegrees(p)).toArray()));
            System.out.println("------------------------------------");
            System.out.printf("Advance ratio J         = %.3f%n", J);
            System.out.printf("Number of blades B        = %d%n", bestB);
            System.out.printf("Propeller radius (m)       = %.2f (Diameter = %.2f m)%n",
                    bestRadius, 2.0 * bestRadius);
            System.out.printf("Chord scale factor         = %.2f%n", bestChordScale);
            System.out.printf("Shaft speed RPM            = %d rpm%n", bestRPM);
            System.out.printf("Tip Mach number            = %.3f%n",
                    (bestRPM * Math.PI / 30.0) * bestRadius / aSound);
            System.out.println("------------------------------------");
            System.out.printf("Total Thrust T             = %.2f N%n", bestThrust);
            System.out.printf("Total Torque Q             = %.2f N·m%n", bestTorque);
            System.out.printf("Required Power             = %.2f W%n", bestPower);
            System.out.printf("Propulsive Efficiency η    = %.3f%n", bestEfficiency);
        } else {
            // If no combination satisfied the thrust target, inform the user.
            System.out.println("No valid design found within the specified search ranges.");
            System.out.println("Consider expanding the search ranges or refining increments.");
        }
    }
}
