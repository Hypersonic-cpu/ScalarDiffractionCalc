(* Imaging script *)

#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
#load "Discrete.fs"
#load "Optics.fs"

open MathNet.Numerics.LinearAlgebra
open LiquidLens
open LiquidLens.Optics
open MathNet.Numerics.Data.Matlab
open Plotly.NET
open Plotly.NET.LayoutObjects

let pi = System.Math.PI
let rad x = x / 180. * pi
(* Test part - Spherical lens *)
let lens = 
    let var = [|0.0 .. 1e-4 .. pi/12.|]
    let rho = var |> Array.map (fun x -> 2.0*sin x)
    let z   = var |> Array.map (fun x -> 2.0*(1.-cos x))
    let diff = var |> Array.map tan
    Discrete.Curve(rho, z, diff, 8.0)

let src = vector [0.; -1.; 0.] 
let trace = 
    [|- rad 5. .. 1e-3 .. rad 15.|]
    |> Array.map (fun x -> vector [0.; sin x; cos x])
    |> Array.map (fun v -> Discrete.Line(src, v))

let ptrXY = 
    trace 
    |> Array.map (Optical.TryRefractRay 1.5 1.0 lens)
    |> Array.filter (fun x -> x.IsSome)
    |> Array.map (fun x -> x.Value)
    |> Array.map (Discrete.Calculate.LineZPlaneIntercept 20.0)
    |> Array.map (fun v -> v.[0], v.[1])

Chart.Scatter(ptrXY, mode=StyleParam.Mode.Markers, Marker=TraceObjects.Marker.init(Size=2, Opacity=0.95))
|> Chart.withXAxis(
    LinearAxis.init(Range=StyleParam.Range.ofMinMax(-1.5, 1.5))
)
|> Chart.withYAxis(
    LinearAxis.init(Range=StyleParam.Range.ofMinMax(-1.5, 1.5))
)
|> Chart.show
