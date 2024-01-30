module Mandelbrot.Graphics.Render
  ( scanline
  , scanline1
  ) where

import Data.Complex (Complex((:+)))

import Mandelbrot.Numerics (Square, Approx)

import Mandelbrot.Graphics.Transform (Transform(..))
import Mandelbrot.Graphics.Compute (Compute, compute, isInterior)

scanline1 :: (RealFloat r, Square r, Approx r) => Int -> Transform r -> r -> Int -> Int -> [Compute r]
scanline1 width transform er maxiters y = go True [0 .. width - 1]
  where
    go bias (x:xs) =
      let px = compute bias er (forward transform (fromIntegral x :+ fromIntegral y, 1)) maxiters
      in  px : go (isInterior px) xs
    go _ _ = []

scanline :: (RealFloat r, Square r, Approx r) => Int -> Int -> Transform r -> r -> Int -> [[Compute r]]
scanline width height transform er maxiters = map (scanline1 width transform er maxiters) [0 .. height - 1]
