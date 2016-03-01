{-| This module defines common data structures for biosequences,
    i.e. data that represents nucleotide or protein sequences.

    Basically, anything resembling or wrapping a sequence should
    implement the 'BioSeq' class (and 'BioSeqQual' if quality information
    is available).

    The data types are mostly wrappers from lazy bytestrings from
    'Data.ByteString.Lazy' and 'Data.ByteString.Lazy.Char8', but most users
    of this module should not need to access the underlying data types directly.
-}

{-# LANGUAGE GeneralizedNewtypeDeriving, DeriveDataTypeable #-}

module Bio.Core.Sequence (
  -- * Data definitions
  Qual (..), Offset (..),
  SeqData (..), SeqLabel (..), QualData (..),

  -- * Class definitions
  BioSeq (..), BioSeqQual (..),

  -- * Helper functions
  toFasta, toFastaQual, toFastQ,

  module Data.Stringable -- stringable doesn't compile with older GHC

  ) where

import qualified Data.ByteString.Lazy.Char8 as LC
import qualified Data.ByteString.Lazy as L
import Data.Int
import Data.Typeable (Typeable)
import Data.Word
import Data.String
import Data.Stringable hiding (length)
import Data.Monoid

-- | Sequence data are lazy bytestrings of ASCII characters.
newtype SeqData  = SeqData { unSD :: LC.ByteString }
                 deriving (Eq,Ord,IsString,Show,Typeable,Stringable)

instance Monoid SeqData where
  mempty = SeqData mempty
  mappend (SeqData s1) (SeqData s2) = SeqData (mappend s1 s2)
  mconcat = SeqData . mconcat . map unSD

-- | Like sequence data, sequence labels are lazy bytestrings of ASCII characters.
newtype SeqLabel = SeqLabel { unSL :: LC.ByteString }
                 deriving (Eq,Ord,IsString,Show,Typeable,Stringable)

instance Monoid SeqLabel where
  mempty = SeqLabel mempty
  mappend (SeqLabel s1) (SeqLabel s2) = let
    (i1:r1) = LC.words s1
    (i2:r2) = LC.words s2
    sid = mconcat [i1,(LC.pack ":"),i2]
    in SeqLabel (LC.unwords ([sid]++r1++[LC.pack ":"]++r2))
  -- mconcat default

-- | A quality value is in the range 0..255.
newtype Qual     = Qual { unQual :: Word8 }
                 deriving (Show,Eq,Ord,Num,Enum,Real,Integral,Typeable)

-- | Quality data are lazy bytestrings of 'Qual's.
newtype QualData = QualData { unQD :: L.ByteString }
                 deriving (Eq,Ord,Show,Typeable,Stringable)

instance Monoid QualData where
  mempty = QualData mempty
  mappend (QualData s1) (QualData s2) = QualData (mappend s1 s2)
  mconcat = QualData . mconcat . map unQD

-- | An 'Offset' is a zero-based index into a sequence
newtype Offset   = Offset { unOff :: Int64 }
                 deriving (Show,Eq,Ord,Num,Enum,Real,Integral,Typeable)

-- | The 'BioSeq' class models sequence data, and any data object that
--   represents a biological sequence should implement it.
class BioSeq s where
  seqid     :: s -> SeqLabel -- ^ Sequence identifier (typically first word of the header)
  seqid = seqlabel
  seqheader :: s -> SeqLabel -- ^ Sequence header (may contain whitespace), by convention the
                             --   first word matches the 'seqid'
  seqheader = seqlabel
  seqdata   :: s -> SeqData  -- ^ Sequence data
  seqlength :: s -> Offset   -- ^ Sequence length

--  slice     :: s -> (Offset,Offset) -> s  -- ^ Cut a slice of a sequence.
--  copy      :: s -> s                     -- ^ Create a copy of a sequence.

  seqlabel  :: s -> SeqLabel -- ^ Deprecated.  Instead, use 'seqid' if you 
                             -- want the unique ID, or 'seqheader' if you 
                             -- want the FASTA style header with ID and comments.
  seqlabel = seqid 

{-# DEPRECATED seqlabel "Warning: 'seqlabel' is deprecated, use 'seqid' or 'seqheader' instead." #-}

-- | Any 'BioSeq' can be formatted as Fasta, 60-char lines.
toFasta :: BioSeq s => s -> LC.ByteString -- any kind of string-like data type?  Use builder?
toFasta s = LC.concat (gt:unSL (seqheader s):nl:wrap (unSD $ seqdata s))
  where wrap x = if LC.null x then [] else let (ln,rest) = LC.splitAt 60 x in ln : nl : wrap rest
        nl = LC.pack "\n"
        gt = LC.pack ">"

-- | The BioSeqQual class extends 'BioSeq' with quality data.  Any correspondig data object
--   should be an instance, this will allow Fasta formatted quality data 'toFastaQual', as
--   well as the combined FastQ format (via 'toFastQ').
class BioSeq sq => BioSeqQual sq where
  seqqual  :: sq -> QualData

-- | Output Fasta-formatted quality data (.qual files), where quality values are output as
--   whitespace-separated integers.
toFastaQual :: BioSeqQual s => s -> LC.ByteString
toFastaQual s = LC.concat (gt:unSL (seqheader s):nl:wrap (L.unpack $ unQD $ seqqual s))
  where wrap x = if null x then [] else let (ln,rest) = splitAt 20 x in LC.pack (unwords $ map show ln) : nl : wrap rest
        nl = LC.pack "\n"
        gt = LC.pack ">"

-- | Output FastQ-formatted data.  For simplicity, only the Sanger quality format is supported,
--   and only four lines per sequence (i.e. no line breaks in sequence or quality data).
toFastQ :: BioSeqQual s => s -> LC.ByteString
toFastQ s = LC.unlines [LC.cons '@' (unSL $ seqid s)
                       , unSD (seqdata s)
                       , LC.cons '+' (unSL $ seqid s)
                       , L.map (+33) (unQD $ seqqual s)]
