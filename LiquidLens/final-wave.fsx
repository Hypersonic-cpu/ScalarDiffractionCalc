#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
;;
#load "Library.fs"
#load "Discrete.fs"
#load "Optics.fs"
#load "Surface.fs"
;;
open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.IntegralTransforms
open MathNet.Numerics.Data.Matlab
open Plotly.NET
open Plotly.NET.LayoutObjects

open LiquidLens.Surface
open LiquidLens.Discrete
open LiquidLens.Optics

let mm = 1e-3
let nm = 1e-9
let um = 1e-6
let pi = System.Math.PI
let J_ = complex 0. 1.
let round (x:float) (n:int) = System.Math.Round(x,n)
let deg rad = rad / pi * 180.
let rad deg = deg / 180. * pi
(* Segment 0. Parameters *)

;;
Control.NativeProviderPath <- @"C:\Users\74674\.nuget\packages\mathnet.numerics.mkl.win-x64\3.0.0\runtimes\win-x64\native\";;
Control.UseNativeMKL()
;;

(* Plot *)
let oResize (mat:Matrix<complex>) =
    let N = mat.ColumnCount
    Array2D.init (N/8-1) (N/8-1) (fun i j -> mat.[8*i+4, 8*j+4])
    |> Matrix.Build.DenseOfArray
let qResize (mat:Matrix<complex>) =
    let N = mat.ColumnCount
    Array2D.init (N/4-1) (N/4-1) (fun i j -> mat.[N/2 + (i-N/8) * 4, N/2 + (j-N/8) * 4])
    |> Matrix.Build.DenseOfArray
let plotPhase (name:string) (mat:Matrix<complex>) =
    Chart.Heatmap(
        mat 
            |> Matrix.map(fun x -> x.Phase) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.withTitle (name+" phase")
    |> Chart.show
let plotMag (name:string) (mat:Matrix<complex>) =
    Chart.Heatmap(
        mat 
            |> Matrix.map(fun x -> x.Magnitude) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.withTitle (name+" mag")
    |> Chart.show
;;

let waveLen = 507. * nm
let waveNum = 2.*pi/waveLen
let nWater = 1.33

let solvePSF radius objDist bias N saveIt = 
    let fileName = $"D:/MatlabDrive/LiqLens/Data/0F-Drop-{radius}.mat"

    let aperture = float radius * um
    let curvatureC:float = 
        let pC = MatlabReader.Read(fileName, "p_c").[0,0]
        pC / Phys.Gamma / 2.
    let focalLen = ( curvatureC*(nWater-1.) )**(-1.)
    let sptlSize = 10.0*mm
    let outputSize = 0.0288*sptlSize
    // let objDist = 20.5*mm
    // 1./objDist + 1./imgDist = 1. / focalLen
    let imgDist = (1./focalLen - 1./objDist) ** (-1.)
    // let bias = 1.6*mm

    let lensThickness = 
        let data = MatlabReader.Read(fileName, "rho_z_theta")
        let lensZMax = -data.[1, data.ColumnCount-1]
        Interpolation
            .CubicSpline
            .InterpolateHermiteSorted(
                data.Row(0).AsArray(), 
                data.Row(1).AsArray(), 
                data.Row(2).AsArray() |> Array.map tan  
            )
            .Interpolate
        |> fun itpl -> (itpl >> (fun x -> x+lensZMax))

    printfn $"Sample point N required: {Optical.WaveLeastSampleNumber waveLen sptlSize}"
    printfn $"U {round (objDist*1000.) 2}\tV {round (imgDist*1000.) 2}\tF {round (focalLen*1000.) 2}"
    // let N = 4096 // stdin.ReadLine() |> int |> fun x -> 2*x/2
    // let freqRange = float N / sptlSize
    // let frequencySample index = (float index / (float N-1.)-0.5) * freqRange
    let spatialSample index = (float index / (float N-1.)-0.5) * sptlSize

    let waveFunc =
        (
        Array2D.init N N
            (fun i j -> 
                let fi, fj = (spatialSample i, spatialSample j)
                Optical.PSFUniform waveLen objDist (fi+bias) fj
                |> Optical.LensPhaseTransform waveNum lensThickness nWater aperture fi fj
            )
        |> Matrix.Build.DenseOfArray
        )
    // let backupWave = Matrix.Build.DenseOfMatrix waveFunc
    // waveFunc |> oResize |> plotPhase "After lens"
    // waveFunc |> oResize |> plotMag "After lens"

    waveFunc |> Optical.TransferWaveInPlace sptlSize waveLen (1.00 * imgDist)
    waveFunc |> Calculate.ResampleInPlace (0., 0.) sptlSize (bias/objDist*imgDist, 0.) outputSize // (0.39*(objDist/imgDist)*sptlSize)
    // waveFunc |> oResize |> plotMag "80% img dist"

    // waveFunc |> Optical.TransferWaveInPlace (0.30*(objDist/imgDist)*sptlSize) waveLen (0.20 * imgDist)
    // waveFunc |> Calculate.ResampleInPlace (bias/objDist*0.80*imgDist, 0.) (0.30*(objDist/imgDist)*sptlSize) (bias/objDist*imgDist, 0.) (0.02496*sptlSize)
    waveFunc |> oResize |> plotMag "100% img dist"

    // let filePath = $"D:/MatlabDrive/LiqLens/Data/0F-PSF-C-r{radius}-bias{round (atan(bias/objDist)*1e2) 0}.mat"
    let filePath = $"D:/MatlabDrive/LiqLens/Data/0F-PSF-V-r{radius}-u{round (objDist*1000.0) 0}.mat"
    printfn $"Save to\n\t{filePath} \n\t" // ? [Y|N]"
    // match stdin.ReadLine() with
    // | "Y" -> 
    if saveIt then MatlabWriter.Write(filePath, [| waveFunc |> oResize; matrix [[complex outputSize 0.]]; matrix [[complex imgDist 0.]]; matrix [[complex objDist 0.]]|], [| "psf"; "edge_len"; "img_dis"; "obj_dis"|] )
    // | _ -> ()
;;

for uDis in [| 26.*mm .. 3.*mm .. 35.*mm|] do
    printfn $"obj dis {uDis * 1e3} mm"
    solvePSF 4900 uDis 0. 6400 true
// solvePSF 4900 (27.4*mm) (0.*um) 6400 true
#quit;;

// let RADIUS = 3200
// let uDis =
//     let fileName = $"D:/MatlabDrive/LiqLens/Data/0F-Drop-{RADIUS}.mat"
//     // let aperture = float RADIUS * um
//     let curvatureC:float = 
//         let pC = MatlabReader.Read(fileName, "p_c").[0,0]
//         pC / Phys.Gamma / 2.
//     2.0 * ( curvatureC*(nWater-1.) )**(-1.)
// let M = 6400
// let bDis = uDis * tan (vector [|0.00..0.02..0.08|])
// for bias in bDis do
//     solvePSF RADIUS uDis bias M true
// #quit
