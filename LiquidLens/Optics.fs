namespace LiquidLens.Optics

open LiquidLens.Discrete
open MathNet.Numerics
open MathNet.Numerics.IntegralTransforms
open MathNet.Numerics.LinearAlgebra

[<AbstractClass; Sealed>]
type Optical =
    /// <summary>
    /// Point spread function of a point source in a uniform space.
    /// </summary>
    /// <param name="waveLength"></param>
    /// <param name="z">Distance from source to observe plane.</param>
    /// <param name="x"></param>
    /// <param name="y"></param>
    static member PSFUniform waveLength (z:float) (x:float) (y:float) = 
        let rho = sqrt (z**2. + x**2. + y**2.)
        let waveNumber = 2.0 * System.Math.PI / waveLength
        exp(complex 0. 1. * waveNumber * rho)
                * z/rho * (1. + 1./complex 0. 1./waveNumber/rho) / complex 0. 1. / waveLength / rho
    
    static member LensPhaseTransform (waveNumber:float) (lensThickness: float->float) nIndex aperture x y (value:complex) = 
        let r = sqrt (x**2. + y**2.)
        if r < aperture then value * exp(complex 0. 1. * waveNumber * lensThickness r * (nIndex-1.))
        else complex 0. 0.

    static member TryRefractRay nRear nFront (rotSurf: Curve) (rayIn: Line) =
        match Calculate.TryRotSurfLineIntercept rotSurf rayIn with
        | Some (intercept) ->
            // Must in front -> rear direction.
            let normalVec = Calculate.RotSurfNormalVector rotSurf intercept
            let angleIn = 
                rayIn.Dir * normalVec 
                |> fun x -> 
                    if x > 1. then 1.
                    elif x < -1. then -1. 
                    else x
                |> acos
            let auxVec = 
                rayIn.Dir - (rayIn.Dir * normalVec) * normalVec
                |> fun x -> 
                    if x.L2Norm()<>0. then x / x.L2Norm()
                    else vector [0.;0.;0.]
            let angleRe = 
                asin (nFront / nRear * sin angleIn)
                |> fun r -> (r+System.Math.PI) % System.Math.PI
            normalVec * cos angleRe + auxVec * sin angleRe
            |> fun dir -> Line(intercept, dir)
            |> fun x -> Some(x)
        | None -> None
    
    /// <summary>
    /// The least sample number to satisfy ``2 * Pi * f_max = k``.
    /// </summary>
    /// <param name="waveLength"> </param>
    /// <param name="spatialFullLen"> Spatial size of input. </param>
    static member WaveLeastSampleNumber waveLength spatialFullLen =
        // Delta f * Delta X = 1 / N.
        // (N * Delta f) * (N * Delta X) = N
        // 2.0 * fMax * spaticalFullLen = N
        // Since fMax = 1./waveLength (k = 2pi*fMax),
        // N = 
        int (spatialFullLen * 2.0 / waveLength)

    // For even N.
    static member private FFTShift (mat: Matrix<'T>) =
        let original = Matrix.Build.DenseOfMatrix mat
        let N =  mat.RowCount
        mat.[N/2.., N/2..] <- original.[..N/2-1, ..N/2-1]   // I
        mat.[..N/2-1, N/2..] <- original.[N/2.., ..N/2-1]   // II
        mat.[..N/2-1, ..N/2-1] <- original.[N/2.., N/2..]   // III
        mat.[N/2.., ..N/2-1] <- original.[..N/2-1, N/2..]   // IV

    static member private WaveTransferFunction (waveNumber: float) (deltaZ:float) fx fy =
        // delta is real.
        let delta = waveNumber**2. - 4.*System.Math.PI**2. * (fx**2.+fy**2.)
        if delta <= 0. then complex 0.0 0.0
        else exp (complex 0. 1. * deltaZ * sqrt(delta))

    /// <summary>
    /// Transfer a wave at z=0 to z=Z.
    /// </summary>
    /// <param name="fullLen"> = y_max - y_min = x_max - x_min in samples. </param>
    /// <param name="waveFunc"> Original wave function of each sampling point. </param>
    /// <returns>
    /// Wave at +Z, a complex matrix.
    /// </returns>
    static member TransferWaveInPlace (fullLen: float) waveLength Z (waveFunc: Matrix<complex>)=
        Control.UseNativeMKL()
        let N = waveFunc.RowCount
        let frequencyRange = float N / fullLen
        let frequencySample index = (float index / (float N-1.)-0.5) * frequencyRange
        let transfer = Optical.WaveTransferFunction (2.*System.Math.PI/waveLength) Z
        // printfn $"{N}\t{frequencyRange}\t{frequencyRange/float N * fullLen}"
        waveFunc |> Fourier.Forward2D
        waveFunc |> Optical.FFTShift
        for fx in 0..N-1 do
            for fy in 0..N-1 do
            waveFunc.[fx, fy] <- 
                waveFunc.[fx, fy] * transfer (frequencySample fx) (frequencySample fy)
        waveFunc |> Optical.FFTShift
        waveFunc |> Fourier.Inverse2D
        ()
