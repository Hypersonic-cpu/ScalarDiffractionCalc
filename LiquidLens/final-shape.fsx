#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
;;
#load "Library.fs"
#load "Surface.fs"
#load "Discrete.fs"
#load "Optics.fs"
;;
open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.IntegralTransforms
open MathNet.Numerics.Data.Matlab
open Plotly.NET
open Plotly.NET.LayoutObjects

open LiquidLens.Discrete
open LiquidLens.Surface
open LiquidLens.Optics

let mm = 1e-3
let nm = 1e-9
let um = 1e-6
let pi = System.Math.PI
let J_ = complex 0. 1.
let round (x:float) (n:int) = System.Math.Round(x,n)
let deg rad = rad / pi * 180.
let rad deg = deg / 180. * pi
;;
(* Segment 0. Parameters *)
// Default at 25.C
// let waterSurfGamma = 74e-3
// let contactAngle = rad 81.
;;
(* Plot *)
let oResize (mat:Matrix<complex>) =
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
;;

let N = 1000
let radiusMsr: float array = [|3.8*mm|] // [|2.4*mm .. 0.4*mm .. 6.0*mm|]
let dropMsr = radiusMsr |> Array.map Drop
let pCenter, matrixReturn= 
    dropMsr
    |> Array.map(fun x -> x.determineFullData 5. 75. (x.Radius*3.) N )
    |> Array.unzip
    |> fun (pC, matArr) ->
        let cutOffIndex = 
            matArr 
            |> Array.mapi(fun i m -> 
                m.Row(2).AsArray() 
                |> Array.findIndex(fun angle -> angle < dropMsr.[i].ContactAngle ) 
            )
        pC,
        (matArr, cutOffIndex)
        ||> Array.map2(fun mat index -> mat.[0.., ..index])

(* Checking Image. *)
(radiusMsr, matrixReturn)
||> Array.iter2(
    fun r (mat:Matrix<float>) ->
        Chart.Scatter(mat.Row(0), mat.Row(1), mode=StyleParam.Mode.Lines)
        |> Chart.withXAxis(
            LinearAxis.init(
                Title = Title.init($"Radius = {round (r*1000.0) 3} mm"), 
                Range = StyleParam.Range.ofMinMax(0., 5.0e-3)
                )
        )
        |> Chart.withYAxis(
            LinearAxis.init(
                Range = StyleParam.Range.ofMinMax(-5.0e-3, 0.)
                )
        )
        |> Chart.show
) 
printfn $"Save? [Y|N]"
match stdin.ReadLine() with 
| "Y" -> 
    radiusMsr
    |> Array.iteri(fun idx r ->
        let fileName = $"D:/MatlabDrive/LiqLens/Data/0F-Drop-{round (1e6*r) 0}.mat"
        printfn $"  [Y|N] ready to write {fileName}"
        if stdin.ReadLine()="Y" then
            MatlabWriter.Write(fileName,
                [| matrix [[r]]; matrix [[pCenter.[idx]]]; matrixReturn.[idx]|],
                [| "r_msr"; "p_c"; "rho_z_theta" |]
            )
    )
| _ -> ()

// #quit;;
