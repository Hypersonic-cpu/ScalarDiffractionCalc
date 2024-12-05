namespace LiquidLens.Surface

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.OdeSolvers
open System.Collections

type R  = float
type VR = Vector<float>

/// <summary>
/// Physics constant in SI unit.
/// </summary>
[<Sealed; AbstractClass>]
type Phys =
    static member Gamma = 74.62e-3 // 301 K.
    static member Density = 996.23 // 301 K.
    static member Gravity = 9.80

type Drop (radius: float) =
    member val Radius = radius
    member val ContactAngle = - System.Math.PI * 81. / 180.
        with get, set
    /// <summary>
    /// Tolerance length (SI meter).
    /// </summary>
    member val LengthTolerance = 2e-8
        with get, set
    /// <summary>
    /// Tolerance angle (rad).
    /// </summary>
    member val AngleTolerance = 1e-6
        with get, set
    
    member val MaximumIterations = 100
        with get, set
    
    /// <summary>
    /// d [rho; z; theta] / ds
    /// </summary>
    /// <param name="pCenter">Pressure difference at rho=0.</param>
    /// <param name="v">[rho; z; theta]</param>
    member this.diff (pCenter:R) (s: R) (v:VR) =
        [|
            cos v.[2];
            sin v.[2];        
            if abs s < this.LengthTolerance
            then -0.5 * pCenter / Phys.Gamma
            else 
            - sin v.[2] / v.[0] 
               - (pCenter - Phys.Density * Phys.Gravity * v.[1]) / Phys.Gamma
        |]
        |> vector

    /// <summary>
    /// Start [rho; z; theta] when DeltaP0 = pCenter.
    /// </summary>
    /// <param name="pCenter"></param>
    /// <param name="rangeS"></param>
    member this.TryStartFromP0 (pCenter, rangeS, N) =
        RungeKutta.FourthOrder(vector [|0.; 0.; 0.|], 
            0., rangeS, N, this.diff pCenter)

    /// <summary>
    /// Determine pressure difference pCenter * maximum of param S.
    /// </summary>
    /// <param name="maxValueS">Possible maximum of param S.</param>
    /// <returns> Pressure difference at center. </returns>
    member this.determineParam pMin pMax maxValueS N =
        let radiusAtPC pCenter = 
            // printfn $"Now pC = {pCenter}"
            let odeRes =
                this.TryStartFromP0(pCenter, maxValueS, N)
                |> Matrix.Build.DenseOfColumnVectors
            Interpolation
                .CubicSpline
                .InterpolateAkima(odeRes.Row(2), odeRes.Row(0))
                .Interpolate(this.ContactAngle)
            
        let objF = radiusAtPC >> fun x -> x - this.Radius
        let mutable pl = pMin
        let mutable pr = pMax
        let mutable maximumIteration = 0
        while abs (objF pl) > this.LengthTolerance && (maximumIteration <= this.MaximumIterations)do
            let mid = (pl + pr) / 2. 
            if objF mid * objF pl < 0. then pr <- mid
            else pl <- mid
            maximumIteration <- maximumIteration+1
            // printfn $"{pl}\t{mid}\t{pr}"
        pl

    /// <summary>
    /// Determine the shape of this drop.
    /// </summary>
    /// <param name="maxValueS">Possible maximum of param S.</param>
    /// <returns> pCenter * curveXY </returns>
    member this.determineShape pCMin pCMax maxParamS N =
        let pCenter = this.determineParam pCMin pCMax maxParamS N 
        pCenter,
        this.TryStartFromP0(pCenter, maxParamS, N)
        |> Array.map(fun v -> v[0], v[1])

    /// <summary>
    /// Determine full shape data of this drop.
    /// </summary>
    /// <param name="pCenterRange">Possible range for pCenter.</param>
    /// <param name="maxValueS">Possible maximum of param S.</param>
    /// <returns> pCenter * Matrix [rho(row0), z(row1), theta(row2)] </returns>
    member this.determineFullData pCMin pCMax maxParamS N =
        let pCenter = this.determineParam pCMin pCMax maxParamS N 
        pCenter,
        this.TryStartFromP0(pCenter, maxParamS, N)
        |> Matrix.Build.DenseOfColumnVectors
