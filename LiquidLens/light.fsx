#I "C:/Users/74674/.nuget/ref/"
#load "MathNet.fsx"
#load "Plotly.NET.fsx"
#load "Discrete.fs"
// #load "Surface.fs"
#load "Optics.fs"

open MathNet.Numerics.LinearAlgebra
open LiquidLens
// open LiquidLens.Surface
open LiquidLens.Optics
open MathNet.Numerics.Data.Matlab
open Plotly.NET

let pi = System.Math.PI

// module GeometryTest =
//     let curveS = vector [0.0 .. 0.001 .. 0.4] * pi
//     let curveX = curveS |> Vector.map sin
//     let curveY = curveS |> Vector.map cos
//     let curveD = curveS |> Vector.map (fun x -> -tan x)
//     let curve = Discrete.Curve (curveX, curveY, curveD, 0.0)
//     let ptArr = [
//         vector [0.;0.];
//         vector [0.144114446;0.2177799999777777];
//         vector [0.144214446;0.2177799999777777];
//         vector [0.144314446;0.2177799999777777];
//         vector [0.144414446;0.2177799999777777];
//         vector [0.144514446;0.2177799999777777];
//         vector [0.144614446;0.2177799999777777];
//     ]
//     // x^2 + y^2 + z^2 = 1.
//     let correct = 
//         ptArr
//         |> List.map(fun x ->
//                 let v = vector [2.*x[0]; 2.*x[1]; 2.* sqrt (1. - x.[0]**2. - x.[1]**2.)]
//                 v / v.L2Norm()
//             )
//     let calc = 
//         ptArr
//         |> List.map 
//             (Discrete.Calculate.RotSurfNormalVector curve)
// ;;
// GeometryTest.correct;;
// GeometryTest.calc;;

module PrecisionTest =
    let GenerateErr ds iota =
        // let ln = Discrete.Line(vector [0.;0.;0.], vector [1.; 0. ;1.])
        printfn $"Gen strt {System.DateTime.Now}"
        let T0 = System.DateTime.Now
        let var = [0.0..ds..pi/12.]
        (*
         *  Test Case: Spherical Surf.
         *  rho = 2 sin x
         *  rho'= 2 cos x
         *  z   = 2 (1 - cos x)
         *  z'  = 2 sin x
         *  dz/d rho = tan x
         *)
        let vRho = 
            var 
            |> List.map (fun x -> 2.0 * sin x) 
            |> vector
        let vZ = 
            var 
            |> List.map (fun x -> 2.0 * (1. - cos x)) 
            |> vector
        let vDiff =
            var
            |> List.map tan 
            |> vector
        printfn $"Gen list ends {(System.DateTime.Now-T0).TotalMilliseconds}ms"
        let surf = Discrete.Curve(
            vRho,
            vZ,
            vDiff, 
            8.0
        )
        printfn $"Gen curv ends {(System.DateTime.Now-T0).TotalMilliseconds}ms"

        let lineArr = 
            let theta = pi/4.0
            Discrete.Line(vector [0.;0.;0.], 
                vector [sin iota * cos theta; sin iota * sin theta; cos iota]) 

        printfn $"Gen line ends {(System.DateTime.Now-T0).TotalMilliseconds}ms"
        let lineArrR =
            lineArr
            |> Optical.TryRefractRay 1.5 1.0 surf 

        printfn $"Refract ends {(System.DateTime.Now-T0).TotalMilliseconds}ms"
        lineArrR
            |> fun x -> x.Value
            |> fun v -> - v.Point.[0] / v.Dir.[0] + v.Point.[2]
            |> fun x -> x-20.

    let dsVec = vector [1e-3; 1e-4; 1e-5; 1e-6]
    let diVec = 1e-6 * vector [6. .. 10.1 .. 1506.]
    let errAtDs ds  =
        diVec 
        |> Vector.map(fun di -> GenerateErr ds di)
        |> fun y -> 
            Chart.Scatter(x=diVec, y=y, mode=StyleParam.Mode.Lines_Markers)
        |> Chart.withTitle $"log={round (log ds/ log 10.)}"
        |> Chart.show
    dsVec   
    |> Vector.iter errAtDs
    // behave well when rho > 50e-6

// module SphShape =
//     let mutable ds = 1e-2
//     printfn $"Gen strt {System.DateTime.Now}"
//     let T0 = System.DateTime.Now
//     let var = [0.0..ds..pi/12.]
//     (*
//      *  Test Case: Spherical Surf.
//      *  rho = 2 sin x
//      *  rho'= 2 cos x
//      *  z   = 2 (1 - cos x)
//      *  z'  = 2 sin x
//      *  dz/d rho = tan x
//      *)
//     let vRho = 
//         var 
//         |> List.map (fun x -> 2.0 * sin x) 
//         |> vector
//     let vZ = 
//         var 
//         |> List.map (fun x -> 2.0 * (1. - cos x)) 
//         |> vector
//     let vDiff =
//         var
//         |> List.map tan 
//         |> vector    
//     let surf = Discrete.Curve(
//         vRho,
//         vZ,
//         vDiff,
//         0.0
//     )
    
// module RefractTest =
//     let var = vector [0.0 .. 0.1 .. 1.0]
//     let rho = var |> Vector.map (fun x -> x)
//     let z = var |> Vector.map (fun x  -> -x)
//     let diff = var |> Vector.map (fun _ -> -1.)
//     // plane: z = 1.0
//     let pln = Discrete.Curve(rho, z, diff, 0.2)
//     let nrm x = Discrete.Calculate.RotSurfNormalVector pln x
//     let itc x = Discrete.Calculate.TryRotSurfLineIntercept pln x 
//     // ray (0,0,0) -> (1,0,1)
//     let ray = Discrete.Line(vector [0.;0.;0.], vector [1.;0.;1.])
//     let ray'= Optical.TryRefractRay 1.5 1.0 pln ray
