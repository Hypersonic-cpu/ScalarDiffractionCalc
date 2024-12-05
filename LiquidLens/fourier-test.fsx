#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
// #load "Discrete.fs"
// #load "Optics.fs"

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics
open MathNet.Numerics.Data.Matlab
open MathNet.Numerics.IntegralTransforms

open Plotly.NET
open Plotly.NET.LayoutObjects
;;
Control.NativeProviderPath <- @"C:\Users\74674\.nuget\packages\mathnet.numerics.mkl.win-x64\3.0.0\runtimes\win-x64\native\";;
Control.UseNativeMKL();;
let pi = System.Math.PI
let J_ = complex 0. 1.
;;

let N = 512
let waveLength = 500e-9     // SI. 500 nm.
let waveNumber = 1. / waveLength
let waveRadius = 3e-3       // SI. 3 mm.

let BiLinearInterpolation (mat:Matrix<float>) (totalRange:float) x y =
    let fM = float mat.RowCount
    let fi, fj = 
        (fM-1.) * (x/totalRange+0.5),
        (fM-1.) * (y/totalRange+0.5)
    let il = floor fi 
    let iu = if ceil fi = il then il+1. else ceil fi
    let jl = floor fj 
    let ju = if ceil fj = jl then jl+1. else ceil fj
    // printf  $"I\t{int il}\t{fi}\t{int iu}\t\t"
    // printfn  $"J\t{int jl}\t{fj}\t{int ju}\t\t"
    // printfn $"{mat.[int il, int jl]}\t{mat.[int il, int ju]}"
    // printfn $"{mat.[int iu, int jl]}\t{mat.[int iu, int ju]}"
    // printfn ""
    (
        (fi-il) * mat.[int iu, int jl] + (iu-fi) * mat.[int il, int jl],
        (fi-il) * mat.[int iu, int ju] + (iu-fi) * mat.[int il, int ju]
    )
    |> fun (lV, uV) -> (fj-jl) * uV + (ju-fj) * lV

let generate waveRange  =
// let waveRange  = 100e-3     // SI. .5mm. Edge length of square sample area.
    (*
    *   /
    *  /
    *  |          .
    *  |          ^ z = R. phase=0 here.
    *  \
    *   \
    *)
    let f x y = 
        let rho = sqrt (waveRadius**2. + x**2. + y**2.)
        exp(J_*waveNumber*rho)
            * waveRadius/rho * (1. + 1./J_/waveNumber/rho) / J_ / waveLength / rho
    let spatialX = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * waveRange|]
    let spatialY = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * waveRange|]
    let waveFunc =
        Array2D.init N N 
            // x = i / (N-1) * W - W/2
            (fun i j ->
                (spatialX.[i], spatialY.[j]) ||> f
            )
        |> Matrix.Build.DenseOfArray
    Chart.Heatmap(
        waveFunc 
            |> Matrix.map(fun x -> x.Magnitude) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.show
    waveFunc

let transfer waveIn waveRange outRange dz =
    let sptlSample = waveRange / (float N - 1.)
    let freqSample = 1./float N/sptlSample          // Delta f * Delta x = 1/N
    let freqRange  = (float N-1.) * freqSample
    let freqncyX = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * freqRange|]
    let freqncyY = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * freqRange|]
    let waveFunc: Matrix<complex> = Matrix.Build.DenseOfMatrix waveIn
    Fourier.Forward2D(waveFunc, FourierOptions.Default)
    let FFTShift (mat:Matrix<'T>) =
        let original = Matrix.Build.DenseOfMatrix mat
        let N = mat.ColumnCount
        mat.[(N+1)/2.., (N+1)/2..] <- original.[..(N+1)/2-1, ..(N+1)/2-1]   // I
        mat.[..(N+1)/2-1, (N+1)/2..] <- original.[(N+1)/2.., ..(N+1)/2-1]   // II
        mat.[..(N+1)/2-1, ..(N+1)/2-1] <- original.[(N+1)/2.., (N+1)/2..]   // III
        mat.[(N+1)/2.., ..(N+1)/2-1] <- original.[..(N+1)/2-1, (N+1)/2..]   // IV
    FFTShift waveFunc

    // Chart.Heatmap(
    //     waveFunc 
    //         |> Matrix.map(fun x -> x.Magnitude) 
    //         |> (fun x -> Matrix.toRowArrays x)
    // )
    // |> Chart.show

    let TransferH (z:float) fx fy =
        // delta is real.
        let delta = waveNumber**2. - 4.*pi**2. * (fx**2.+fy**2.)
        if delta <= 0. then complex 0.0 0.0
        else exp (J_*z*sqrt(delta))

    let waveTransfer = TransferH (dz)
    
    let waveFunc' =
        waveFunc
        |> Matrix.mapi(fun fx fy value -> 
            value * 
                waveTransfer (freqncyX.[fx]) (freqncyY.[fy]))
    FFTShift waveFunc'
    
    Fourier.Inverse2D(waveFunc', FourierOptions.Default)
    // Chart.Heatmap(
    //     waveFunc'
    //         |> Matrix.map(fun x -> x.Magnitude) 
    //         |> (fun x -> Matrix.toRowArrays x)
    // )
    // |> Chart.show

    let sampleL i = (float i / (float N-1.)-0.5) * outRange
    let magnitudeIntpl = BiLinearInterpolation (waveFunc' |> Matrix.map(fun x -> x.Magnitude))  waveRange
    let realIntpl = BiLinearInterpolation (waveFunc' |> Matrix.map(fun x -> x.Real))            waveRange
    let imagIntpl = BiLinearInterpolation (waveFunc' |> Matrix.map(fun x -> x.Imaginary))       waveRange
    let interpolateIt i j = 
        let fi, fj = sampleL i, sampleL j
        // try
        let aux = complex (realIntpl fi fj) (imagIntpl fi fj)
        let mag = magnitudeIntpl fi fj
        complex (mag * cos aux.Phase) (mag * sin aux.Phase)
        // with _ ->
        //     printfn $"{i} {j}\t\t{fi}\t{fj}\t{waveRange}"
        //     raise(System.Exception("______________________-"))
    Array2D.init N N 
        (fun i j -> interpolateIt i j)
    |> Matrix.Build.DenseOfArray

let showIt (matr:Matrix<complex>) = 
    Chart.Heatmap(
        matr
            |> Matrix.map(fun x -> x.Magnitude) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.show

let mutable ranges = 2e-3
let m0 = generate ranges
let m1 = transfer m0 ranges 1.8e-3 (1.0*waveRadius)
showIt m1
// let m2 = transfer m1 150e-3 100e-3 (0.2*waveRadius)
// showIt m2
// let m3 = transfer m2 100e-3 50e-3 (0.2*waveRadius)
// showIt m3
// let m4 = transfer m3 50e-3 20e-3 (0.2*waveRadius)
// showIt m4
// let m5 = transfer m4 20e-3 10e-3 (0.2*waveRadius)
// showIt m5

