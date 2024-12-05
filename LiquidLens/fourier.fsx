#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
#load "Discrete.fs"
#load "Optics.fs"

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.IntegralTransforms
open MathNet.Numerics.Data.Matlab
open Plotly.NET
open Plotly.NET.LayoutObjects

open LiquidLens.Optics
open LiquidLens.Discrete

Control.NativeProviderPath <- @"C:\Users\74674\.nuget\packages\mathnet.numerics.mkl.win-x64\3.0.0\runtimes\win-x64\native\";;
Control.UseNativeMKL()
let pi = System.Math.PI
let J_ = complex 0. 1.
;;

let mm = 1e-3
let nm = 1e-9
let waveLength = 595.*nm
let waveNumber = 2.*pi/waveLength
// Delta f * N = 1. / Delat X = const.
let sptlRange = 6.0*mm
let N = 
    printfn $"{Optical.WaveLeastSampleNumber waveLength sptlRange}"
    stdin.ReadLine() |> int |> fun x -> (x*2+1)/2
let freqRange = float N / sptlRange
let frequencySample index = (float index / (float N-1.)-0.5) * freqRange
let spatialSample index = (float index / (float N-1.)-0.5) * sptlRange

(* Plot *)
let qResize (mat:Matrix<complex>) =
    let N = mat.ColumnCount
    Array2D.init (N/8-1) (N/8-1) (fun i j -> mat.[N/2 + (i-N/16) * 8, N/2 + (j-N/16) * 8])
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

(* Generate original wave *)
// Plane - Convex lens.
let focalLen = 9.*mm
let nLen = 1.500
let apertureRadius = 4.0*mm
// 1/f_E = Power = 1/R * (nLen-1)
// f_E * (nLen-1) = R
let PassLens x y value = 
    let radius = focalLen * (nLen-1.)
    let z0 = sqrt (radius**2. - apertureRadius**2.)
    match sqrt (radius**2. - x**2. - y**2.) - z0 with
    | thickness when thickness>0. ->
        value * exp(complex 0. 1. * waveNumber * thickness * (nLen-1.))
    | _ -> complex 0. 0.
let PassLenIndex i j value =
    let x = spatialSample i 
    let y = spatialSample j
    PassLens x y value

let BiasX = 1.5*mm
let waveFunc = 
    Array2D.init N N (fun i j -> 
        Optical.PSFUniform waveLength (2.0*focalLen) 
            (spatialSample i + BiasX) (spatialSample j) )
    |> Array2D.mapi PassLenIndex
    |> Matrix.Build.DenseOfArray
waveFunc |> qResize |> plotMag      "Rear lens surface -"
//waveFunc |> qResize |> plotPhase    "Rear lens surface -"

Optical.TransferWaveInPlace sptlRange waveLength (1.20*focalLen) waveFunc
//waveFunc |> qResize |> plotMag      "Interstage -"
;;
Calculate.ResampleInPlace (0., 0.) sptlRange (0.2*BiasX, 0.) (0.60*sptlRange) waveFunc
//waveFunc |> qResize |> plotMag      "Interstage -"
;;
Optical.TransferWaveInPlace (0.60*sptlRange) waveLength (0.60*focalLen) waveFunc
//waveFunc |> qResize |> plotMag      "Interstage2-"
;;
Calculate.ResampleInPlace (0.2*BiasX, 0.) (0.60*sptlRange) (0.75*BiasX, 0.) (0.2*sptlRange) waveFunc
waveFunc |> qResize |> plotMag      "Interstage2-"
;;
Optical.TransferWaveInPlace (0.2*sptlRange) waveLength (0.20*focalLen) waveFunc
waveFunc |> qResize |> plotMag      "Imaging plane -"
;;
//waveFunc |> qResize |> plotPhase    "Imaging plane -"

#quit;;
