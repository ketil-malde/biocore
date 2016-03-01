{- | Utility functions, tables, and data structures for manipulating sequences. -}
module Bio.Core.Sequence.Util (
  -- * Functions for nucleotide sequences
  revcompl, revqual, compl
  -- * Functionality for amino acid sequences                     
  , iupac, iupac', toIUPAC, fromIUPAC, Amino(..)
  -- * Translation of sequences
  , trans_tbl, translate  
  -- * Amino acid properties             
  , Property(..), property_list             
  ) where

import Bio.Core.Sequence
import qualified Data.ByteString.Lazy.Char8 as LC
import Data.List (unfoldr)
import Data.Char (toUpper)
import Data.Maybe (fromJust)

-- | Calculate the reverse-complement of a sequence. This only makes
--   sense for nucleotide data, of course.
revcompl :: SeqData -> SeqData
revcompl = SeqData . LC.reverse . LC.map compl . unSD

-- | Reverse QualData (to be used in conjunction with 'revcompl')
revqual :: QualData -> QualData
revqual = QualData . LC.reverse . unQD

-- | The complement of a nucleotide, including IUPAC wildcards.
compl :: Char -> Char
compl 'A' = 'T'
compl 'T' = 'A'
compl 'C' = 'G'
compl 'G' = 'C'
compl 'a' = 't'
compl 't' = 'a'
compl 'c' = 'g'
compl 'g' = 'c'
compl 'n' = 'n'
compl 'N' = 'N'
compl 'R' = 'Y' --	A or G
compl 'Y' = 'R' --	C or T
compl 'S' = 'S' --	G or C
compl 'W' = 'W' --	A or T
compl 'K' = 'M' --	G or T
compl 'M' = 'K' --	A or C
compl 'B' = 'V' --	C or G or T
compl 'D' = 'H' --	A or G or T
compl 'H' = 'D' --	A or C or T
compl 'V' = 'B' --	A or C or G
compl x = x     --      including '-' for gaps

-- -- | Read the character at the specified position in the sequence.
-- {-# INLINE (!) #-}
-- (!) :: Sequence -> Offset -> Char
-- (!) (Seq _ (SeqData bs) _) = LC.index bs . unOff

------------------------------------------------------------
-- Amino acid (protein) stuff
------------------------------------------------------------

-- | Translate a nucleotide sequence into the corresponding protein
--   sequence.  This works rather blindly, with no attempt to identify ORFs
--   or otherwise QA the result.
translate :: SeqData -> Offset -> [Amino]
translate (SeqData s') o' = unfoldr triples (s',o')
   where triples (s,o) = 
             if unOff o > LC.length s - 3 then Nothing
             else Just (trans1 (map (LC.index s . unOff) [o,o+1,o+2]),(s,o+3))

-- prop_trans s = tail (translate s 0) == translate s 3

trans1 :: String -> Amino
trans1 = maybe Xaa id . flip lookup trans_tbl . map (repUT . toUpper)
    where repUT x = if x == 'U' then 'T' else x -- RNA uses U for T

-- | List of the standard amino acids.
data Amino = Ala | Arg | Asn | Asp | Cys | Gln | Glu | Gly
           | His | Ile | Leu | Lys | Met | Phe | Pro | Ser
           | Thr | Tyr | Trp | Val 
           | STP | Asx | Glx | Xle | Xaa -- unknowns
     deriving (Show,Eq)

-- | The standard codon table for amino acid translation.  Others exist, notably for
--   mitochondria and some prokaryotes or archea.
trans_tbl :: [(String,Amino)]
trans_tbl = [("GCT",Ala),("GCC",Ala),("GCA",Ala),("GCG",Ala),
             ("CGT",Arg),("CGC",Arg),("CGA",Arg),("CGG",Arg),("AGA",Arg),("AGG",Arg),
             ("AAT",Asn),("AAC",Asn),
             ("GAT",Asp),("GAC",Asp),
--             ("RAT",Asx),("RAC",Asx), -- IUPAC: R is purine (A or G)
             ("TGT",Cys),("TGC",Cys),
             ("CAA",Gln),("CAG",Gln),
             ("GAA",Glu),("GAG",Glu),
--             ("SAA",Glx),("SAG",Glx), -- IUPAC: S is C or G
             ("GGT",Gly),("GGC",Gly),("GGA",Gly),("GGG",Gly),
             ("CAT",His),("CAC",His),
             ("ATT",Ile),("ATC",Ile),("ATA",Ile),
             ("TTA",Leu),("TTG",Leu),("CTT",Leu),("CTC",Leu),("CTA",Leu),("CTG",Leu),
             ("AAA",Lys),("AAG",Lys),
             ("ATG",Met),
             ("TTT",Phe),("TTC",Phe),
             ("CCT",Pro),("CCC",Pro),("CCA",Pro),("CCG",Pro),
             ("TCT",Ser),("TCC",Ser),("TCA",Ser),("TCG",Ser),("AGT",Ser),("AGC",Ser),
             ("ACT",Thr),("ACC",Thr),("ACA",Thr),("ACG",Thr),
             ("TAT",Tyr),("TAC",Tyr),
             ("TGG",Trp),
             ("GTT",Val),("GTC",Val),("GTA",Val),("GTG",Val),
             ("TAG",STP),("TGA",STP),("TAA",STP)
            ]
-- todo: extend with more IUPAC nucleotide wildcards?

-- | Convert a list of amino acids to a sequence in IUPAC format.
toIUPAC :: [Amino] -> SeqData
toIUPAC = SeqData . LC.pack . map (fromJust . flip lookup iupac)

-- | Convert a sequence in IUPAC format to a list of amino acids.
fromIUPAC :: SeqData -> [Amino]
fromIUPAC = map (maybe Xaa id . flip lookup iupac' . toUpper) . LC.unpack . unSD

-- | Table mapping amino acids to their IUPAC letters.
iupac :: [(Amino,Char)]
iupac = [(Ala,'A')        ,(Arg,'R')        ,(Asn,'N')
        ,(Asp,'D')        ,(Cys,'C')        ,(Gln,'Q')
        ,(Glu,'E')        ,(Gly,'G')        ,(His,'H')
        ,(Ile,'I')        ,(Leu,'L')        ,(Lys,'K')
        ,(Met,'M')        ,(Phe,'F')        ,(Pro,'P')
        ,(Ser,'S')        ,(Thr,'T')        ,(Tyr,'Y')
        ,(Trp,'W')        ,(Val,'V')
        ,(Asx,'B') -- Asn or Asp
        ,(Glx,'Z') -- Gln or Glu
        ,(Xle,'J') -- Ile or Leu
        ,(STP,'X')        ,(Xaa,'X')
        ]

-- | Table mapping IUPAC characters to amino acids (the inverse of 'iupac').
iupac' :: [(Char,Amino)]
iupac' = map (\(a,b)->(b,a)) iupac

data Property = Small | Polar | Hydrophobic 
              | Aliphatic | Aromatic
              | Tiny | Charged | Positive | Negative
              deriving (Eq,Ord,Read,Show,Enum)

property_list :: [(Property, String)]
property_list = 
  [ (Small,"AGCSPNDTVC")
  , (Tiny,"AGCS")
  , (Negative,"DE")
  , (Positive,"KRH")
  , (Charged,"DEKRH")
  , (Polar,"NSCTYWQDEKRH")
  , (Aromatic,"FYHW")
  , (Aliphatic,"ILV")
  , (Hydrophobic,"FYHWILVCGATMK")
  ]
