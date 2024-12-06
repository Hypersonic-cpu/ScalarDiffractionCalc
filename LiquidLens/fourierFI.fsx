#load "../LibLoader.fsx"
#load "../MklLoader.fsx"

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.IntegralTransforms
open Plotly.NET
open Plotly.NET.LayoutObjects

Control.NativeProviderPath <- MklLoader.MKLProviderPath
Control.UseNativeMKL()
let pi = System.Math.PI
let J_ = complex 0. 1.
;;

(* Settings *)
let waveLength = 0.0008     // mm. = 800 nm. (NIR)
let waveNumber = 2.*pi/waveLength
// Delta f * N = 1. / Delat X = const.
let freqRange = 2.0 / waveLength
let sptlRange = 2.
let N = int (freqRange * sptlRange)

let waveRadius = 0.05          // mm. 3 mm.

printfn $"N {N}"

stdin.ReadLine() |> ignore

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

let FFTShift (mat:Matrix<'T>) =
    let original = Matrix.Build.DenseOfMatrix mat
    let N = mat.ColumnCount
    mat.[(N+1)/2.., (N+1)/2..] <- original.[..(N+1)/2-1, ..(N+1)/2-1]   // I
    mat.[..(N+1)/2-1, (N+1)/2..] <- original.[(N+1)/2.., ..(N+1)/2-1]   // II
    mat.[..(N+1)/2-1, ..(N+1)/2-1] <- original.[(N+1)/2.., (N+1)/2..]   // III
    mat.[(N+1)/2.., ..(N+1)/2-1] <- original.[..(N+1)/2-1, (N+1)/2..]   // IV

(* Generate original wave *)
let g x y = 
    let rho = sqrt (waveRadius**2. + x**2. + y**2.)
    exp(J_*waveNumber*rho)
        * (-waveRadius)/rho * (1. + 1./J_/waveNumber/rho) / J_ / waveLength / rho
// let g x y = 
//     let rho = sqrt (waveRadius**2. + x**2. + y**2.)
//     exp(J_*waveNumber*rho) * (waveRadius)/rho / J_ / waveLength / rho
// let f x y =
//     if x**2. + y**2. < 5e-6 then complex 1e6 0.
//     // if abs x < 0.1e-3 && abs y < 2e-3 then complex 1e6 0.
//     else complex 0. 0.
let spatialX = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * sptlRange|]
let spatialY = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * sptlRange|]
let freqncyX = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * freqRange|]
let freqncyY = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * freqRange|]
// Spherical wave at z=radius from point impulse.
let waveFunc =
    Array2D.init N N 
        // x = i / (N-1) * W - W/2
        (fun i j ->
            (spatialX.[i], spatialY.[j]) ||> g
        )
    |> Matrix.Build.DenseOfArray

// let waveFunc2 =
//     Array2D.init N N 
//         // x = i / (N-1) * W - W/2
//         (fun i j ->
//             (spatialX.[i], spatialY.[j]) ||> f
//         )
//     |> Matrix.Build.DenseOfArray

waveFunc |> qResize |> plotPhase "Z0"
waveFunc |> qResize |> plotMag "Z0"
// stdin.ReadLine() |> ignore

(* Forward FFT *)
Fourier.Forward2D(waveFunc, FourierOptions.Default)
FFTShift waveFunc
waveFunc |> qResize |> plotPhase "FT shifted"
waveFunc |> qResize |> plotMag "FT shifted"



(* Freq domain transfer *)
let TransferH (z:float) fx fy =
    // delta is real.
    let delta = waveNumber**2. - 4.*pi**2. * (fx**2.+fy**2.)
    if delta <= 0. then complex 0.0 0.0
    else exp (J_*z*sqrt(delta))
let waveTransfer = TransferH waveRadius
let waveTransferI = TransferH (-waveRadius)
let waveFunc' =
    waveFunc
    |> Matrix.mapi(fun fx fy value -> 
        value * 
            waveTransferI (freqncyX.[fx]) (freqncyY.[fy]))
    // |> Matrix.mapi(fun fx fy value -> 
    //     value * 
    //         waveTransferI (freqncyX.[fx]) (freqncyY.[fy]))
(* In INVERSE direction, backtracking -> src. *)
waveFunc' |> qResize |> plotPhase "Transfered"
waveFunc' |> qResize |> plotMag "Transfered"
FFTShift waveFunc'

(* Inverse FFT *)
Fourier.Inverse2D(waveFunc', FourierOptions.Default)
(* SHOULD BE A POINT *)
waveFunc' |> qResize |> plotPhase "Final"
waveFunc' |> qResize |> plotMag "Final"

// Fourier.Forward2D(waveFunc', FourierOptions.Default)
// FFTShift waveFunc'
// let waveFunc'' =
//     waveFunc'
//     |> Matrix.mapi(fun fx fy value -> 
//         value * 
//             waveTransfer (freqncyX.[fx]) (freqncyY.[fy]))
// FFTShift waveFunc''
// Fourier.Inverse2D(waveFunc'', FourierOptions.Default)
// waveFunc'' |> qResize |> plotPhase "Revive"
// waveFunc'' |> qResize |> plotMag "Revive"
