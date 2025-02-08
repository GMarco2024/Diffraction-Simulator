//
//  ContentView.swift
//  Diffraction Simulator
//


import SwiftUI
import Foundation

// MARK: - Input Validation

struct ValidationError: Error {
    let message: String
}

struct SimulationParameters {
    let w0: Double      // initial beam width
    let r0: Double      // initial radius of wavefront curvature
    let ell0: Double    // initial transverse coherence width
    let lambda: Double  // wavelength
    let d: Double       // period of grating
    let nu1: Double     // grating 1 open fraction
    let Z1: Double      // position (mm)
    let theta: Double   // twist between gratings, in degrees
    
    static func validate(
        beamWidth: String,
        curvature: String,
        coherence: String,
        wavelength: String,
        gratingPeriod: String,
        gratingNu1: String,
        gratingZ1: String,
        twist: String
    ) throws -> SimulationParameters {
        guard let w0 = Double(beamWidth), w0 > 0 else {
            throw ValidationError(message: "Beam width must be a positive number")
        }
        guard let r0 = Double(curvature) else {
            throw ValidationError(message: "Curvature must be a valid number")
        }
        guard let ell0 = Double(coherence), ell0 > 0 else {
            throw ValidationError(message: "Coherence must be a positive number")
        }
        guard let lambda = Double(wavelength), lambda > 0 else {
            throw ValidationError(message: "Wavelength must be a positive number")
        }
        guard let d = Double(gratingPeriod), d > 0 else {
            throw ValidationError(message: "Grating period must be a positive number")
        }
        guard let nu1 = Double(gratingNu1), nu1 > 0, nu1 <= 1 else {
            throw ValidationError(message: "Nu1 must be between 0 and 1")
        }
        guard let Z1 = Double(gratingZ1), Z1 > 0 else {
            throw ValidationError(message: "Z1 must be a positive number")
        }
        guard let theta = Double(twist) else {
            throw ValidationError(message: "Twist must be a valid number")
        }
        
        return SimulationParameters(w0: w0, r0: r0, ell0: ell0, lambda: lambda, d: d, nu1: nu1, Z1: Z1, theta: theta)
    }
}

// MARK: - Utility Functions

/// Normalized sinc function.
func sinc(_ x: Double) -> Double {
    return abs(x) < 1e-12 ? 1.0 : sin(x) / x
}

/// Computes the effective propagation distance.
func zp(z: Double, v: Double) -> Double {
    return (v * z) / (v + z)
}

/// Returns a common scaling factor used by wz and ellz.
/// The ratio abs(z/zp) and the square-root term are common for the beam spread calculations.
func scalingFactor(z: Double, r0: Double, ell0: Double, w0: Double, lambda: Double) -> Double {
    let zpVal = zp(z: z, v: r0)
    let ratio = abs(z / zpVal)
    let sqrtTerm = sqrt(1 + pow(lambda * zpVal, 2) / pow(ell0 * w0, 2))
    return ratio * sqrtTerm
}

/// Beam width at propagation distance z.
func wz(z: Double, r0: Double, ell0: Double, w0: Double, lambda: Double) -> Double {
    return w0 * scalingFactor(z: z, r0: r0, ell0: ell0, w0: w0, lambda: lambda)
}

/// Coherence width at propagation distance z.
func ellz(z: Double, r0: Double, ell0: Double, w0: Double, lambda: Double) -> Double {
    return ell0 * scalingFactor(z: z, r0: r0, ell0: ell0, w0: w0, lambda: lambda)
}

/// Computes the radius of curvature at propagation distance z.
func rz(z: Double, r0: Double, ell0: Double, w0: Double, lambda: Double) -> Double {
    let zpVal = zp(z: z, v: r0)
    let factor = 1 + pow(lambda * zpVal, 2) / pow(ell0 * w0, 2)
    let denom = 1.0 - (zpVal / (z * factor))
    return z / denom
}

// MARK: - Core Simulation Functions

/// Computes the intensity profile before the grating.
/// Uses a Gaussian profile shaped by the beam width.
func gp0(z: Double, r0: Double, ell0: Double, w0: Double, lambda: Double,
         xpoints: Int, Xmin: Double, Xmax: Double) -> [(x: Double, value: Double)] {
    let w = wz(z: z, r0: r0, ell0: ell0, w0: w0, lambda: lambda)
    var result = [(x: Double, value: Double)]()
    let dx = (Xmax - Xmin) / Double(xpoints - 1)
    
    for i in 0..<xpoints {
        let x = Xmin + Double(i) * dx
        let value = exp(-Double.pi * pow(x / w, 2))
        result.append((x: x, value: value))
    }
    return result
}

/// Computes the intensity profile after the first grating.
func gp1(z12: Double, r1: Double, ell1: Double, w1: Double, lambda: Double,
         xpoints: Int, Xmin: Double, Xmax: Double, nu1: Double, d: Double) -> [(x: Double, value: Double)] {
    let ell2 = ellz(z: z12, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let w2 = wz(z: z12, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let r2 = rz(z: z12, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let cutoff = 1e-3
    let lim = 4
    let dx = (Xmax - Xmin) / Double(xpoints - 1)
    let xs: [Double] = (0..<xpoints).map { Xmin + Double($0) * dx }
    var intensity = Array(repeating: 0.0, count: xpoints)
    
    for n in (-lim...lim) {
        for m in (-lim...lim) {
            let dn = Double(n - m)
            let dm = (Double(n) + Double(m)) / 2.0
            let baseCoef = sinc(Double(n) * .pi * nu1) * sinc(Double(m) * .pi * nu1) * pow(nu1, 2)
            let envCoef = exp(-Double.pi * pow(dn * lambda * z12 / (d * ell2), 2))
            let coef = baseCoef * envCoef
            if coef >= cutoff {
                for (i, x) in xs.enumerated() {
                    let arg = (x - dm * lambda * z12 / d) / w2
                    let phase = 2 * Double.pi * (dn / d) * (x - dm * lambda * z12 / d) * (1 - z12 / r2)
                    intensity[i] += coef * exp(-Double.pi * pow(arg, 2)) * cos(phase)
                }
            }
        }
    }
    return Array(zip(xs, intensity))
}

/// Computes the intensity profile after the second grating.
func gp2(z12: Double, z23: Double, theta: Double,
         r1: Double, ell1: Double, w1: Double, lambda: Double,
         xpoints: Int, Xmin: Double, Xmax: Double, d: Double) -> [(x: Double, value: Double)] {
    let z13 = z12 + z23
    let ell3 = ellz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let w3 = wz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    // Precompute phase constant that is independent of n and m.
    let phi = 2 * Double.pi * lambda * z23 / (d * d)
    let cutoff = 1e-3
    let lim = 4
    let dx = (Xmax - Xmin) / Double(xpoints - 1)
    let xs: [Double] = (0..<xpoints).map { Xmin + Double($0) * dx }
    var intensity = Array(repeating: 0.0, count: xpoints)
    
    // Convert theta to radians once.
    let thetaRad = theta * Double.pi / 180.0
    
    for n in (-lim...lim) {
        for m in (-lim...lim) {
            let dn = Double(n - m)
            let dm = (Double(n) + Double(m)) / 2.0
            let coef1 = exp(-Double.pi * pow(dn * sin(thetaRad) * lambda * z23 / (d * ell3), 2))
            let coef2 = exp(-Double.pi * pow(lambda * z23 * (dn * cos(thetaRad) + dm * (z13 / z23)) / (d * ell3), 2))
            let coef = coef1 * coef2
            if coef >= cutoff {
                for (i, x) in xs.enumerated() {
                    let arg = (x - (lambda * z23 / d) * (Double(n) * cos(thetaRad) + Double(m) * (z13 / z23))) / w3
                    intensity[i] += coef * exp(-Double.pi * pow(arg, 2)) * cos(phi)
                }
            }
        }
    }
    return Array(zip(xs, intensity))
}

/// Runs the simulation and returns a 2D intensity matrix.
func simulateIntensityMap(params: SimulationParameters) -> [[Double]] {
    // Characteristic length along the propagation axis.
    let LT = pow(params.d, 2) / params.lambda
    let NZ = floor(20e-3 / LT)
    let Z1 = params.Z1
    let Z2 = Z1 + NZ * LT
    let Z3 = Z2 + (Z2 - Z1)
    
    let Zmin: Double = 1e-3
    let Zmax: Double = Z3 + 0.5 * Z1
    let zpoints = 500
    let dZ = (Zmax - Zmin) / Double(zpoints - 1)
    
    let xpoints = 500
    let Xmin = -100 * params.d
    let Xmax = 100 * params.d
    
    // Precompute beam parameters at the grating position.
    let r1 = rz(z: Z1, r0: params.r0, ell0: params.ell0, w0: params.w0, lambda: params.lambda)
    let ell1 = ellz(z: Z1, r0: params.r0, ell0: params.ell0, w0: params.w0, lambda: params.lambda)
    let w1 = wz(z: Z1, r0: params.r0, ell0: params.ell0, w0: params.w0, lambda: params.lambda)
    
    var intensityMatrix = Array(repeating: Array(repeating: 0.0, count: xpoints), count: zpoints)
    
    for i in 0..<zpoints {
        let zpos = Zmin + Double(i) * dZ
        var profile: [(x: Double, value: Double)] = []
        
        if zpos < Z1 {
            profile = gp0(z: zpos, r0: params.r0, ell0: params.ell0, w0: params.w0,
                          lambda: params.lambda, xpoints: xpoints, Xmin: Xmin, Xmax: Xmax)
        } else if zpos < Z2 {
            let z12 = zpos - Z1
            profile = gp1(z12: z12, r1: r1, ell1: ell1, w1: w1,
                          lambda: params.lambda, xpoints: xpoints, Xmin: Xmin, Xmax: Xmax,
                          nu1: params.nu1, d: params.d)
        } else {
            let z12 = Z2 - Z1
            let z23 = zpos - Z2
            profile = gp2(z12: z12, z23: z23, theta: params.theta,
                          r1: r1, ell1: ell1, w1: w1, lambda: params.lambda,
                          xpoints: xpoints, Xmin: Xmin, Xmax: Xmax, d: params.d)
        }
        
        for j in 0..<xpoints {
            intensityMatrix[i][j] = profile[j].value
        }
    }
    return intensityMatrix
}

// MARK: - Views

struct IntensityPlotView: View {
    let matrix: [[Double]]
    
    /// Computes a grayscale color for an intensity value normalized between min and max.
    func colorForIntensity(norm: Double) -> Color {
        return Color(white: norm)
    }
    
    var body: some View {
        GeometryReader { geo in
            let rows = matrix.count
            if rows == 0 {
                EmptyView()
            } else {
                let cols = matrix[0].count
                let cellWidth = geo.size.width / CGFloat(cols)
                let cellHeight = geo.size.height / CGFloat(rows)
                let flat = matrix.flatMap { $0 }
                let maxVal = flat.max() ?? 1.0
                let minVal = flat.min() ?? 0.0
                
                Canvas { context, _ in
                    for i in 0..<rows {
                        for j in 0..<cols {
                            let norm = (matrix[i][j] - minVal) / (maxVal - minVal + 1e-12)
                            let rect = CGRect(x: CGFloat(j) * cellWidth,
                                              y: CGFloat(i) * cellHeight,
                                              width: cellWidth,
                                              height: cellHeight)
                            context.fill(Path(rect), with: .color(colorForIntensity(norm: norm)))
                        }
                    }
                }
                // Rotate so that the diffraction pattern displays correctly.
                .rotationEffect(.degrees(-90))
            }
        }
        .aspectRatio(0.75, contentMode: .fit)
    }
}

struct ContentView: View {
    // Input fields with default values.
    @State private var beamWidth = "5e-06"
    @State private var curvature = "-9.99e+20"
    @State private var coherence = "5e-06"
    @State private var wavelength = "3e-09"
    @State private var gratingPeriod = "1e-06"
    @State private var gratingNu1 = "0.5"
    @State private var gratingZ1 = "0.005"
    @State private var twist = "0.0"
    
    @State private var intensityMatrix: [[Double]] = []
    @State private var errorMessage: String = ""
    @State private var showError = false
    
    
    //Beige color. You can customize whatever color you want
    private let customColor = Color(red: 0.96, green: 0.93, blue: 0.85)
    
    var body: some View {
        VStack(spacing: 20) {
            Text("Diffraction Simulator")
                .font(.title)
                .foregroundColor(customColor)
                .padding()
            
            HStack(spacing: 20) {
                // Group beam parameters.
                parameterGroup(title: "Beam parameters", fields: [
                    ("w0 =", $beamWidth),
                    ("r0 =", $curvature),
                    ("ell0 =", $coherence),
                    ("λ =", $wavelength)
                    
                ])
                
                // Group grating parameters.
                parameterGroup(title: "Grating parameters", fields: [
                    ("d =", $gratingPeriod),
                    ("ν₁ =", $gratingNu1),
                    ("Z₁ =", $gratingZ1),
                    ("θ =", $twist)
                ])
            }
            .padding()
            
            Button(action: runSimulation) {
                Text("Run Simulation")
                    .font(.headline)
                    .foregroundColor(.white)
                    .cornerRadius(5)
            }
            .padding(.horizontal)
            
            if !intensityMatrix.isEmpty {
                IntensityPlotView(matrix: intensityMatrix)
                    .frame(height: 300)
                    .padding()
            }
        }
        .padding()
        .alert("Input Error", isPresented: $showError) {
            Button("OK", role: .cancel) { }
        } message: {
            Text(errorMessage)
        }
    }
    
    /// Groups parameters into a VStack with a title and text fields.
    func parameterGroup(title: String, fields: [(String, Binding<String>)]) -> some View {
        VStack(alignment: .leading, spacing: 10) {
            Text(title)
                .font(.headline)
            ForEach(0..<fields.count, id: \.self) { index in
                HStack {
                    Text(fields[index].0)
                    TextField("", text: fields[index].1)
                        .textFieldStyle(RoundedBorderTextFieldStyle())
                        .frame(width: 120)
                }
            }
        }
    }
    
    /// Runs the simulation; validates input, computes the intensity matrix, and handles errors.
    func runSimulation() {
        do {
            let params = try SimulationParameters.validate(
                beamWidth: beamWidth,
                curvature: curvature,
                coherence: coherence,
                wavelength: wavelength,
                gratingPeriod: gratingPeriod,
                gratingNu1: gratingNu1,
                gratingZ1: gratingZ1,
                twist: twist
            )
            intensityMatrix = simulateIntensityMap(params: params)
        } catch let error as ValidationError {
            errorMessage = error.message
            showError = true
            intensityMatrix = []
        } catch {
            errorMessage = "An unexpected error occurred"
            showError = true
            intensityMatrix = []
        }
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
