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
    let nu2: Double     // grating 2 open fraction
    let X2: Double      // transverse offset of second grating
    let Z1: Double      // position of first grating
    let theta: Double   // twist between gratings, in degrees
    let imageEnabled: Bool // whether to use image-charge Fourier coefficients
    
    static func validate(
        beamWidth: String,
        curvature: String,
        coherence: String,
        wavelength: String,
        gratingPeriod: String,
        gratingNu1: String,
        gratingNu2: String,
        gratingX2: String,
        gratingZ1: String,
        twist: String,
        imageEnabled: Bool
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
            throw ValidationError(message: "Grating 1 fraction must be between 0 and 1")
        }
        guard let nu2 = Double(gratingNu2), nu2 > 0, nu2 <= 1 else {
            throw ValidationError(message: "Grating 2 fraction must be between 0 and 1")
        }
        guard let X2 = Double(gratingX2) else {
            throw ValidationError(message: "Grating 2 offset must be a valid number")
        }
        guard let Z1 = Double(gratingZ1), Z1 > 0 else {
            throw ValidationError(message: "Distance to Grating 1 must be a positive number")
        }
        guard let theta = Double(twist) else {
            throw ValidationError(message: "Twist must be a valid number")
        }
        
        return SimulationParameters(
            w0: w0,
            r0: r0,
            ell0: ell0,
            lambda: lambda,
            d: d,
            nu1: nu1,
            nu2: nu2,
            X2: X2,
            Z1: Z1,
            theta: theta,
            imageEnabled: imageEnabled
        )
    }
}

// MARK: - Utility Functions

/// Normalized sinc function matching Mathematica's Sinc[x] = sin(pi x)/(pi x).
func sinc(_ x: Double) -> Double {
    guard abs(x) >= 1e-12 else { return 1.0 }
    let scaled = Double.pi * x
    return sin(scaled) / scaled
}

struct Complex {
    var re: Double
    var im: Double

    init(re: Double = 0.0, im: Double = 0.0) {
        self.re = re
        self.im = im
    }

    static func *(lhs: Complex, rhs: Complex) -> Complex {
        return Complex(
            re: lhs.re * rhs.re - lhs.im * rhs.im,
            im: lhs.re * rhs.im + lhs.im * rhs.re
        )
    }

    static func /(lhs: Complex, rhs: Double) -> Complex {
        return Complex(re: lhs.re / rhs, im: lhs.im / rhs)
    }

    func conjugate() -> Complex {
        return Complex(re: re, im: -im)
    }

    func maxComponentMagnitude() -> Double {
        return max(abs(re), abs(im))
    }
}

struct FourierCoefficients {
    let coefficients: [Int: Complex]

    func value(for index: Int) -> Complex {
        return coefficients[index] ?? Complex()
    }

    static func compute(nu: Double, d: Double, res: Int = 1000, wedgeAngle: Double = 0, tilt: Double = 0, thick: Double = 0) -> FourierCoefficients {
        let width = nu * d
        let alpha = wedgeAngle * Double.pi / 180
        let beta = tilt * Double.pi / 180

        let step = width / Double(res)
        let minposX: Double
        let maxposX: Double

        if beta >= 0 {
            if beta <= alpha {
                minposX = -(width * cos(beta)) / 2.0 + step
                maxposX = (width * cos(beta)) / 2.0 - step
            } else {
                minposX = -(width * cos(beta)) / 2.0 + step - thick * (tan(alpha) - tan(beta))
                maxposX = (width * cos(beta)) / 2.0 - step + thick * (tan(alpha) - tan(beta))
            }
        } else {
            if abs(beta) <= alpha {
                minposX = -(width * cos(beta)) / 2.0 + step
            } else {
                minposX = -(width * cos(beta)) / 2.0 + step - thick * (tan(alpha) - tan(beta))
            }
            maxposX = (width * cos(beta)) / 2.0 - step
        }

        var coefficients = [Int: Complex]()
        for n in -20...20 {
            coefficients[n] = Complex()
        }

        var position = minposX
        while position <= maxposX - step + 1e-15 {
            for n in -20...20 {
                let fc = 2 * Double.pi * Double(n) * position / d
                if var c = coefficients[n] {
                    c.re += cos(fc)
                    c.im += sin(fc)
                    coefficients[n] = c
                }
            }
            position += step
        }

        for n in -20...20 {
            coefficients[n] = (coefficients[n] ?? Complex()) / Double(res)
        }

        return FourierCoefficients(coefficients: coefficients)
    }
}

func normalizedProfile(_ profile: [(x: Double, value: Double)]) -> [(x: Double, value: Double)] {
    // Sanitize values: replace NaN/Inf with 0 and clamp negatives to 0
    let sanitized = profile.map { pair -> (x: Double, value: Double) in
        let v = pair.value
        if v.isNaN || !v.isFinite { return (x: pair.x, value: 0.0) }
        return (x: pair.x, value: max(0.0, v))
    }
    guard let maxValue = sanitized.map({ $0.value }).max(), maxValue > 0 else {
        return sanitized
    }
    return sanitized.map { (x: $0.x, value: $0.value / maxValue) }
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
         xpoints: Int, Xmin: Double, Xmax: Double, nu1: Double, d: Double,
         imageEnabled: Bool, coefficients: FourierCoefficients?) -> [(x: Double, value: Double)] {
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
            let coefficient: Double
            if imageEnabled, let coeffN = coefficients?.value(for: n), let coeffM = coefficients?.value(for: m) {
                coefficient = (coeffN * coeffM.conjugate()).re
            } else {
                coefficient = sinc(Double(n) * nu1) * sinc(Double(m) * nu1) * pow(nu1, 2)
            }
            let envCoef = exp(-Double.pi * pow(dn * lambda * z12 / (d * ell2), 2))
            let coef = coefficient * envCoef
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
         xpoints: Int, Xmin: Double, Xmax: Double, d: Double,
         nu1: Double, nu2: Double, X2: Double,
         imageEnabled: Bool, coefficients1: FourierCoefficients?, coefficients2: FourierCoefficients?) -> [(x: Double, value: Double)] {
    let z13 = z12 + z23
    let ell3x = ellz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let w3x = wz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let r3x = rz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let ell3y = ellz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let r3y = rz(z: z13, r0: r1, ell0: ell1, w0: w1, lambda: lambda)
    let cutoff = 1e-3
    let lim = 4
    let dx = (Xmax - Xmin) / Double(xpoints - 1)
    let xs: [Double] = (0..<xpoints).map { Xmin + Double($0) * dx }
    var intensity = Array(repeating: 0.0, count: xpoints)
    let thetaRad = theta * Double.pi / 180.0
    let cosTheta = cos(thetaRad)
    let sinTheta = sin(thetaRad)
    let d1 = d
    let d2 = d
    let commonPrefactor = 2 * Double.pi * lambda * z23 / (d1 * d1)

    for m1 in (-lim...lim) {
        for m2 in (-lim...lim) {
            for n1 in (-lim...lim) {
                for n2 in (-lim...lim) {
                    let dn = Double(n1 - n2)
                    let n = Double(n1 + n2) / 2.0
                    let dm = Double(m1 - m2)
                    let m = Double(m1 + m2) / 2.0

                    let coefficient: Complex
                    if imageEnabled,
                       let coeffM1 = coefficients2?.value(for: m1),
                       let coeffM2 = coefficients2?.value(for: m2),
                       let coeffN1 = coefficients1?.value(for: n1),
                       let coeffN2 = coefficients1?.value(for: n2) {
                        let prodM = coeffM1 * coeffM2.conjugate()
                        let prodN = coeffN1 * coeffN2.conjugate()
                        coefficient = prodM * prodN
                    } else {
                        let coefM = sinc(Double(m1) * nu1) * sinc(Double(m2) * nu1)
                        let coefN = sinc(Double(n1) * nu2) * sinc(Double(n2) * nu2)
                        coefficient = Complex(re: coefM * coefN, im: 0.0)
                    }

                    let env1 = exp(-Double.pi * pow(dn * sinTheta * lambda * z23 / (d2 * ell3y), 2))
                    let env2 = exp(-Double.pi * pow(lambda * z23 * (dn * cosTheta + dm * (z13 / z23)) / (d1 * ell3x), 2))
                    let complexCoef = Complex(re: coefficient.re * env1 * env2,
                                              im: coefficient.im * env1 * env2)

                    if abs(complexCoef.re) >= cutoff || abs(complexCoef.im) >= cutoff {
                        let phiConst = dn * n * (1 - z23 / r3x) * cosTheta * cosTheta
                            + dn * n * (1 - z23 / r3y) * sinTheta * sinTheta
                            + dn * m * (1 - z13 / r3x) * cosTheta
                            + dm * n * (1 - z13 / r3x) * cosTheta
                            + dm * m * (z13 / z23) * (1 - z13 / r3x)
                        let phiBase = phiConst * commonPrefactor - 2 * Double.pi * dn * X2 / d2
                        let linearTerm = 2 * Double.pi / d2 * (dn * cosTheta * (1 - z23 / r3x) + dm * (1 - z13 / r3x))
                        let shift = lambda * z23 / d1 * (Double(n) * cosTheta + Double(m) * (z13 / z23))

                        for (i, x) in xs.enumerated() {
                            let phase = phiBase - linearTerm * x
                            let envelope = exp(-Double.pi * pow((x - shift) / w3x, 2))
                            let term = complexCoef.re * cos(phase) - complexCoef.im * sin(phase)
                            intensity[i] += term * envelope
                        }
                    }
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
    let coefficients1 = params.imageEnabled ? FourierCoefficients.compute(nu: params.nu1, d: params.d) : nil
    let coefficients2 = params.imageEnabled ? FourierCoefficients.compute(nu: params.nu2, d: params.d) : nil
    
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
                          nu1: params.nu1, d: params.d,
                          imageEnabled: params.imageEnabled, coefficients: coefficients1)
        } else {
            let z12 = Z2 - Z1
            let z23 = zpos - Z2
            profile = gp2(z12: z12, z23: z23, theta: params.theta,
                          r1: r1, ell1: ell1, w1: w1, lambda: params.lambda,
                          xpoints: xpoints, Xmin: Xmin, Xmax: Xmax,
                          d: params.d, nu1: params.nu1, nu2: params.nu2, X2: params.X2,
                          imageEnabled: params.imageEnabled,
                          coefficients1: coefficients1,
                          coefficients2: coefficients2)
        }
        
        let normalized = normalizedProfile(profile)
        for j in 0..<xpoints {
            intensityMatrix[i][j] = normalized[j].value
        }
    }
    return intensityMatrix
}

// MARK: - Views

struct IntensityPlotView: View {
    let matrix: [[Double]]
    
    /// Computes a grayscale color for an intensity value normalized between min and max.
    func colorForIntensity(norm: Double) -> Color {
        let v = norm.isFinite ? norm : 0.0
        let clamped = min(max(v, 0.0), 1.0)
        return Color(white: clamped)
    }
    
    var body: some View {
        VStack(spacing: 12) {
            // Plot Title
            Text("Simulated Intensity Distribution")
                .font(.headline)
                .fontWeight(.bold)
                .foregroundColor(.primary)
                .padding(.top, 4)
            
            HStack(spacing: 8) {
                // Y-Axis Label
                Text("Transverse Position (X)")
                    .font(.subheadline)
                    .foregroundColor(.secondary)
                    .rotationEffect(.degrees(-90))
                    // .fixedSize() ensures the text frame doesn't get clipped after rotation
                    .fixedSize()
                
                // The Plot
                GeometryReader { geo in
                    let rows = matrix.count
                    if rows == 0 {
                        EmptyView()
                    } else {
                        let cols = matrix[0].count
                        
                        let cellWidth = geo.size.width / CGFloat(rows)
                        let cellHeight = geo.size.height / CGFloat(cols)
                        
                        let flat = matrix.flatMap { $0 }
                        let maxVal = flat.max() ?? 1.0
                        let minVal = flat.min() ?? 0.0
                        
                        Canvas { context, _ in
                            for i in 0..<rows {
                                for j in 0..<cols {
                                    let denom = max(maxVal - minVal, 1e-12)
                                    let raw = (matrix[i][j] - minVal) / denom
                                    let norm = raw.isFinite ? raw : 0.0
                                    
                                    let rect = CGRect(
                                        x: CGFloat(i) * cellWidth,
                                        y: CGFloat(cols - 1 - j) * cellHeight,
                                        width: cellWidth,
                                        height: cellHeight
                                    )
                                    context.fill(Path(rect), with: .color(colorForIntensity(norm: norm)))
                                }
                            }
                        }
                    }
                }
                .border(Color.gray.opacity(0.3), width: 1) // Add a border so the edges are visible
            }
            
            // X-Axis Label
            Text("Propagation Distance (Z)")
                .font(.subheadline)
                .foregroundColor(.secondary)
                .padding(.bottom, 4)
        }
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
    @State private var gratingNu2 = "0.5"
    @State private var gratingX2 = "0.0"
    @State private var gratingZ1 = "0.005"
    @State private var twist = "0.0"
    @State private var imageEnabled = false
    @State private var isRunning = false
    @State private var runTask: Task<Void, Never>? = nil
    
    @State private var intensityMatrix: [[Double]] = []
    @State private var errorMessage: String = ""
    @State private var showError = false
    
    // Beige color. You can customize whatever color you want.
    private let customColor = Color(red: 0.96, green: 0.93, blue: 0.85)
    
    var body: some View {
        HStack(spacing: 0) {
            // Left panel: inputs and controls
            VStack(alignment: .leading, spacing: 16) {
                Text("Wave Interference & Diffraction Simulator")
                    .font(.title2)
                    .fontWeight(.bold)
                    .foregroundColor(customColor)
                    .padding(.bottom, 8)

                ScrollView {
                    VStack(spacing: 20) {
                        // Parameter groups in a vertical stack
                        parameterGroup(title: "Initial Beam Parameters", fields: [
                            ("Initial Beam Width (w₀)", $beamWidth),
                            ("Wavefront Curvature Radius (r₀)", $curvature),
                            ("Transverse Coherence Width (ℓ₀)", $coherence),
                            ("Wavelength (λ)", $wavelength)
                        ])

                        parameterGroup(title: "Grating Parameters", fields: [
                            ("Grating Period / Spacing (d)", $gratingPeriod),
                            ("Grating 1 Open Fraction (ν₁)", $gratingNu1),
                            ("Grating 2 Open Fraction (ν₂)", $gratingNu2),
                            ("Grating 2 Offset (X₂)", $gratingX2),
                            ("Distance to Grating 1 (Z₁)", $gratingZ1),
                            ("Twist Angle Between Gratings (θ°)", $twist)
                        ])

                        Toggle("Use image-charge Fourier coefficients", isOn: $imageEnabled)
                            .toggleStyle(SwitchToggleStyle(tint: .blue))
                            .foregroundColor(customColor)
                            .padding(.top, 4)
                    }
                }

                Button(action: runSimulation) {
                    Text(isRunning ? "Running..." : "Run Simulation")
                        .font(.headline)
                        .foregroundColor(.white)
                        .frame(maxWidth: .infinity)
                        .padding(.vertical, 12)
                        .background(isRunning ? Color.gray : Color.blue)
                        .cornerRadius(10)
                }
                .disabled(isRunning)
                .padding(.top, 8)
            }
            .padding()
            .frame(minWidth: 320, idealWidth: 380, maxWidth: 420, maxHeight: .infinity, alignment: .topLeading)
            .background(Color.black.opacity(0.02))

            Divider()

            // Right panel: visualization
            VStack(spacing: 12) {
                if isRunning {
                    ProgressView("Computing diffraction pattern...")
                        .padding(.top)
                }

                if !intensityMatrix.isEmpty {
                    IntensityPlotView(matrix: intensityMatrix)
                        .padding()
                        .frame(maxWidth: .infinity, maxHeight: .infinity)
                } else {
                    // Placeholder when no data
                    ZStack {
                        Color.clear
                        Text("Run the simulation to see results")
                            .foregroundColor(.secondary)
                    }
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                }
            }
            .padding()
            .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .top)
        }
        .padding()
        .alert("Input Error", isPresented: $showError) {
            Button("OK", role: .cancel) { }
        } message: {
            Text(errorMessage)
                .foregroundColor(customColor)
        }
    }
    
    /// Groups parameters into a nicely styled VStack with a title and text fields.
    func parameterGroup(title: String, fields: [(String, Binding<String>)]) -> some View {
        VStack(alignment: .leading, spacing: 14) {
            Text(title)
                .font(.headline)
                .foregroundColor(customColor)
                .padding(.bottom, 2)
            
            ForEach(Array(fields.enumerated()), id: \.0) { idx, field in
                HStack {
                    Text(field.0)
                        .foregroundColor(customColor)
                        .font(.subheadline)
                        .frame(maxWidth: .infinity, alignment: .leading)
                    
                    TextField("", text: field.1)
                        .textFieldStyle(RoundedBorderTextFieldStyle())
                        .frame(width: 110)
                        #if os(iOS)
                        .keyboardType(.numbersAndPunctuation)
                        #endif
                }
            }
        }
        .padding()
        .background(Color.black.opacity(0.04)) // Gives a slight card-like appearance
        .cornerRadius(12)
    }
    
    /// Runs the simulation; validates input, computes the intensity matrix, and handles errors.
    func runSimulation() {
        // Prevent overlapping runs
        if isRunning {
            return
        }

        do {
            let params = try SimulationParameters.validate(
                beamWidth: beamWidth,
                curvature: curvature,
                coherence: coherence,
                wavelength: wavelength,
                gratingPeriod: gratingPeriod,
                gratingNu1: gratingNu1,
                gratingNu2: gratingNu2,
                gratingX2: gratingX2,
                gratingZ1: gratingZ1,
                twist: twist,
                imageEnabled: imageEnabled
            )
            // Cancel any previous pending task just in case
            runTask?.cancel()

            intensityMatrix = []
            isRunning = true
            showError = false

            runTask = Task.detached(priority: .userInitiated) {
                let matrix = simulateIntensityMap(params: params)
                if Task.isCancelled { return }
                await MainActor.run {
                    intensityMatrix = matrix
                    isRunning = false
                }
            }
        } catch let error as ValidationError {
            errorMessage = error.message
            showError = true
            intensityMatrix = []
            isRunning = false
        } catch {
            errorMessage = "An unexpected error occurred"
            showError = true
            intensityMatrix = []
            isRunning = false
        }
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
