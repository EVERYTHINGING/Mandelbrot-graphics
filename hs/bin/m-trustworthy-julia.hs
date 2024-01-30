module Main (main) where

import Control.Monad (forM_)
import Data.Foldable (toList)
import Data.Maybe (mapMaybe)
import System.Environment (getArgs)
import qualified Data.Set as S

import Mandelbrot.Graphics.Trustworthy.Julia
import Mandelbrot.Graphics.Trustworthy.Quadratic
import qualified Mandelbrot.Graphics.BoundingBox.D1 as D1
import Mandelbrot.Graphics.BoundingBox.D2
import Mandelbrot.Graphics.QuadTree

main :: IO ()
main = do
  args <- getArgs
  case args of
    [outfile, sre, sim, sdepth] -> do
      re <- readIO sre
      im <- readIO sim
      depth <- readIO sdepth
      let r :: Double
          r = max (sqrt $ re^2 + im^2) 2
          images = julia (quadratic (singleton re im)) r
          content = svg r (images !! depth)
      writeFile outfile content

svg :: Real r => r -> Image r -> String
{-# SPECIALIZE svg :: Double -> Image Double -> String #-}
{-# SPECIALIZE svg :: Rational -> Image Rational -> String #-}
svg r i = unlines $
  [ "<?xml version='1.0' encoding='UTF-8' ?>"
  , "<svg xmlns='http://www.w3.org/2000/svg' version='1.0' viewBox='" ++ f (-r) ++ " " ++ f (-r) ++ " " ++ f (2 * r) ++ " " ++ f (2 * r) ++ "' width='4096' height='4096'>"
  , "<g stroke='none'>"
  ] ++ map (uncurry s) (toList i) ++
  [ "</g>"
  , "</svg>"
  ]
  where
    f z = show (realToFrac z :: Double)
    s b c = svgBox b (rgb c)
    rgb Black  = "black"
    rgb Gray   = "red"
    rgb White  = "white"
    rgb Marked = "yellow"

svgBox :: Real r => Box r -> String -> String
{-# SPECIALIZE svgBox :: Box Double -> String -> String #-}
{-# SPECIALIZE svgBox :: Box Rational -> String -> String #-}
svgBox (Box (D1.Box xlo xhi) (D1.Box ylo yhi)) fill = unwords $
      [ "<path d='M"
      , x0, ",", y0, "L"
      , x0, ",", y1, "L"
      , x1, ",", y1, "L"
      , x1, ",", y0, "L"
      , x0, ",", y0, "Z' fill='" ++ fill ++ "' />"
      ]
      where
        [x0, y0, x1, y1] = map f [xlo, ylo, xhi, yhi]
        f z = show (realToFrac z :: Double)
