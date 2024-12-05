namespace LiquidLens.Discrete
open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra

type R  = float
type VR = Vector<float>

type Line =
    val Point: VR 
    val Dir:   VR
    new (point: VR, dir: VR) = { Point = point; Dir = dir/dir.L2Norm() }

type Curve =
    val Zmin:   R
    val Zmax:   R
    val OriginZ:R
    /// <summary>
    /// Create new curve for rotational surface.
    /// </summary>
    /// <param name="rho"> Increasinglly sorted requied. </param>
    /// <param name="diff"> dz/d(rho) </param>
    new (rho: VR, z: VR, diff: VR, originZ: R) = 
        {
            Zmin = min z.[0] z.[z.Count-1];
            Zmax = max z.[0] z.[z.Count-1];
            OriginZ = originZ;
            DerivativeAtP =
                Interpolation
                    .CubicSpline
                    .InterpolateBoundariesSorted(
                            rho.AsArray(), diff.AsArray(), 
                            Interpolation.SplineBoundaryCondition.FirstDerivative, 
                            0.0,
                            Interpolation.SplineBoundaryCondition.ParabolicallyTerminated,
                            rho.[rho.Count-1]
                    )
                    .Interpolate;
            InterpolateAtP = 
                Interpolation
                    .CubicSpline
                    .InterpolateHermiteSorted(
                        rho.AsArray(), 
                        z.AsArray(), 
                        diff.AsArray()
                    )
                    .Interpolate
        }
    new (rho: R seq, z: R seq, diff: R seq, originZ: R) = 
        Curve(vector rho, vector z, vector diff, originZ)
    /// <summary> Interpolate z value at rho. </summary>
    val InterpolateAtP  : R -> R
    /// <summary> dz / d rho at rho. </summary>
    val DerivativeAtP   : R -> R

// Math service module.
[<AbstractClass; Sealed>]
type Calculate =
    static member Eps = 0.5e-9

    /// <summary>
    /// Intercept of a rotational surface with a line.
    /// Return None if intercept is out of range.
    /// </summary>
    static member TryRotSurfLineIntercept (rotSurf: Curve) (line:Line) =
        let paramMax, paramMin =
            (rotSurf.Zmax + rotSurf.OriginZ - line.Point.[2]) / line.Dir.[2] + Calculate.Eps,
            (rotSurf.Zmin + rotSurf.OriginZ - line.Point.[2]) / line.Dir.[2] - Calculate.Eps
        // printfn $"{paramMax}\t{paramMin}"
        let objDeltaZ param =
            line.Point + param * line.Dir
            |> fun x -> x 
            |> fun (v:VR) -> sqrt (v.[0]**2. + v.[1]**2.), v.[2]
            |> fun (rhoLine, zLine) 
                -> rotSurf.InterpolateAtP rhoLine + rotSurf.OriginZ - zLine
                    // printfn $"rho {rhoLine}, zLine {zLine}, var {var}"; var
        match
            FindRoots.ofFunction paramMin paramMax objDeltaZ
        with
        | Some par -> Some (line.Point + line.Dir * par)
        | None -> None

    static member BiLinearInterpolation (mat:Matrix<complex>) (totalRange:float) x y =
        let fM = float mat.RowCount
        let fi, fj = 
            (fM-1.) * (x/totalRange+0.5),
            (fM-1.) * (y/totalRange+0.5)
        let il = floor fi 
        let iu = if ceil fi = il then il+1. else ceil fi
        let jl = floor fj 
        let ju = if ceil fj = jl then jl+1. else ceil fj
        let interpolate (func:complex -> float) =
            (
                ((fi-il) * mat.[int iu, int jl] |> func) + ((iu-fi) * mat.[int il, int jl] |> func),
                ((fi-il) * mat.[int iu, int ju] |> func) + ((iu-fi) * mat.[int il, int ju] |> func)
            )
            |> fun (lV, uV) -> (fj-jl) * uV + (ju-fj) * lV
        let magIt = interpolate (fun x -> x.Magnitude)
        let phaseIt = 
            (
                complex 
                    (interpolate (fun x -> x.Real))
                    (interpolate (fun x -> x.Imaginary))
            ).Phase
        complex (magIt * cos phaseIt) (magIt * sin phaseIt)

    /// <summary>
    /// Intercept of a line with plane z=zVal
    /// </summary>
    static member LineZPlaneIntercept (zVal:R) (line:Line) =
        (zVal - line.Point.[2]) / line.Dir.[2]
        |> fun param -> line.Point + param * line.Dir

    /// <summary> Normal vector for rotational surface. </summary>
    /// <param name="rotSurf"> generatrix of surface </param>
    /// <param name="point"> R^2 or R^3 vector (ignore Z-dim) </param>
    static member RotSurfNormalVector (rotSurf: Curve) (point: VR) = 
        let rhoPt = (point.[0]**2. + point.[1]**2.) |> sqrt
        if rhoPt=0. then [ 0.;0.;1. ]
        else
            [
                - rotSurf.DerivativeAtP rhoPt * point.[0] / rhoPt;
                - rotSurf.DerivativeAtP rhoPt * point.[1] / rhoPt;
                1.
            ]
        |> vector
        |> fun x -> x / x.L2Norm()
    
    static member ResampleInPlace inCenter inLen outCenter outLen (matx:Matrix<complex>) =
        let N = matx.RowCount
        let inX, inY = inCenter
        let outX, outY = outCenter
        let sampleL i = (float i / (float N-1.)-0.5) * outLen
        let intpl = Calculate.BiLinearInterpolation matx inLen
        Array2D.init N N 
            (fun i j -> intpl (sampleL i + outX - inX) (sampleL j + outY - inY))
        |> Matrix.Build.DenseOfArray
        |> Matrix.iteri (fun i j v -> matx.[i, j] <- v)
