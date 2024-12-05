#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
#load "Discrete.fs"
#load "Surface.fs"

open MathNet.Numerics.LinearAlgebra
open LiquidLens.Surface
open Plotly.NET
open Plotly.NET.Interactive
open Plotly.NET.LayoutObjects
open MathNet.Numerics.Data.Matlab

module DisplayShape = 
    let r = 0.002410288065766
    let drp = Drop(r)

    let mutable pC = 180.
    while pC < 181. do
        drp.TryStartFromP0(pC, 0.045, 7500)
        |> Array.map (fun x -> x.[0], x.[1])
        |> Array.unzip
        |> fun (rho, z) -> Chart2D.Chart.Scatter(rho, z, mode=StyleParam.Mode.Markers)
        |> Chart.withXAxis(
            LinearAxis.init(
                Title = Title.init($"{System.Math.Round(pC, 3)}"), 
                Range = StyleParam.Range.ofMinMax(0., 5.5e-3)
                )
            )
        |> Chart.withYAxis(
            LinearAxis.init(
                Range = StyleParam.Range.ofMinMax(-5.5e-3, 0.0)
                )
            )
        |> Chart.show
        pC <- pC + 10.

module DetermineShape =
    let r = 0.003_311
    let drp = Drop(r)

    let pCenter, fullData = 
        drp.determineFullData 0. 60. 10e-3 2001
    Chart.Scatter(fullData.Row(0), fullData.Row(1), mode=StyleParam.Mode.Markers)
    |> Chart.show
    let valPCenter = matrix [[pCenter]]

    printfn $"Save? [Y|N]"
    match stdin.ReadLine().[0] with 
    | 'Y' ->
        MatlabWriter.Write(@"D:\MatlabDrive\LiqLens\Data\ex_drop_r33_rebuild.mat", 
                            [|
                                valPCenter; 
                                fullData.[0..0, 0..];
                                fullData.[1..1, 0..];
                                fullData.[2..2, 0..];
                            |], 
                            [|"p_center"; "curve_rho"; "curve_z"; "curve_theta"|]
        )
    | _ -> ()

// module FocalLenWrtRadius =
//     let radii = [|0.1e-3..0.05e-3..4.5e-3|]
//     let pCenter = 
//         radii
//         |> Array.map(fun r ->
//             let drp = Drop(r)
//             drp.determineParam 0.0 3000.0 (r * 3.16) 4501
//             // Chart.Scatter(curveXY, mode=StyleParam.Mode.Markers)
//             // |> Chart.withXAxis(LinearAxis.init(Title=Title.init($"R = {r}\tP = {pC}")))
//             // |> Chart.show
//         )
//     (*Modify Required.*)
//     let height = 
//         Array.map3(fun r pc s -> Drop(r).TryStartFromP0(pc, s, 4501)) radii pCenter sParam
//         |> Array.map(fun va -> va.[va.Length-1].[1])
//     let dataMat = Matrix.Build.DenseOfRowArrays [radii; pCenter; height]
//     Chart.Scatter(radii, height, mode=StyleParam.Mode.Markers)
//     |> Chart.show
//     stdin.ReadLine() |> ignore
//     MatlabWriter.Write(@"D:\MatlabDrive\LiqLens\Data\p_center_vs_radii.mat", dataMat, "radii_p0")

// module CmpWithSphere =
//     let r = 0.002_410_288065766
//     let drp = Drop(r)
//     let pCenter, curveXY, maxParamS = drp.DetermineShape((0., 60.), 10e-3, 10001)
//     let rCurvCircle = 1. / (pCenter / Phys.Gamma / 2.)
//     let sphereXY = 
//         [|
//             for theta in -System.Math.PI/4. .. 0.001 .. 0.
//                 ->  - 2. * rCurvCircle * sin theta * cos theta, 
//                     - 2. * rCurvCircle * sin theta * sin theta
//         |]
//     [
//         Chart.Scatter(curveXY, mode=StyleParam.Mode.Lines)
//         Chart.Scatter(sphereXY, mode=StyleParam.Mode.Lines)
//     ]
//     |> Chart.combine
//     |> Chart.show
//     let matCurve = 
//         curveXY
//         |> Array.unzip
//         |> fun (rho, z) -> Matrix.Build.DenseOfRowArrays [|rho; z|]
//     let valPCenter = matrix [[pCenter]]
//     let valSParam = matrix [[maxParamS]]
//     printfn $"Save? [Y|N]"
//     match stdin.ReadLine().[0] with 
//     | 'Y' ->
//         MatlabWriter.Write(@"D:\MatlabDrive\LiqLens\Data\ex_drop_r33_theory.mat", 
//                             [|valPCenter; matCurve; valSParam|], [|"p_center"; "curve_rho_z"; "max_s_param"|])
//     | _ -> ()
