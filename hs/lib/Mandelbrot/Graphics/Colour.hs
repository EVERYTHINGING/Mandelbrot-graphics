module Mandelbrot.Graphics.Colour
  ( monochrome
  , grayscale
  ) where

import Data.Word (Word8)

import Mandelbrot.Graphics.Compute (Compute(..))

monochrome :: RealFloat r => Compute r -> Word8
monochrome (Exterior{ computeDE = de }) = round (255 * tanh de)
monochrome (Interior{ computeDE = de }) = round (255 * tanh de)
monochrome _ = 0

grayscale :: Word8 -> (Word8, Word8, Word8)
grayscale g = (g, g, g)
