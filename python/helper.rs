use std::*;
use std::collections::HashMap;

"
Helper Classes for Configuration Interaction methods

References:
- Equations from [Szabo:1996]
";
const __authors__: _ = "Tianyuan Zhang";
let __credits__ = vec!["Tianyuan Zhang", "Jeffrey B. Schriber", "Daniel G. A. Smith"];
const __copyright__: _ = "(c) 2014-2018, The Psi4NumPy Developers";
const __license__: _ = "BSD-3-Clause";
const __date__: _ = "2017-05-26";
use itertools::{combinations};
struct Determinant {
alphaObtBits: ST0,
betaObtBits: ST1,
}

impl Determinant {
"
    A class for a bit-Determinant.
    ";
fn __init__<T0, T1, T2, T3>(&self, alphaObtBits: T0, betaObtBits: T1, alphaObtList: T2, betaObtList: T3)  {
"
        Constructor for the Determinant
        ";
if alphaObtBits == 0&&alphaObtList != None {
alphaObtBits = Determinant::obtIndexList2ObtBits(alphaObtList);
}
if betaObtBits == 0&&betaObtList != None {
betaObtBits = Determinant::obtIndexList2ObtBits(betaObtList);
}
self.alphaObtBits = alphaObtBits;
self.betaObtBits = betaObtBits;
}
fn getNumOrbitals<RT>(&self) -> RT {
"
        Return the number of orbitals (alpha, beta) in this determinant
        ";
return (Determinant::countNumOrbitalsInBits(self.alphaObtBits), Determinant::countNumOrbitalsInBits(self.betaObtBits));
}
fn getOrbitalIndexLists<RT>(&self) -> RT {
"
        Return lists of orbital index
        ";
return (Determinant::obtBits2ObtIndexList(self.alphaObtBits), Determinant::obtBits2ObtIndexList(self.betaObtBits));
}
fn getOrbitalMixedIndexList<RT>(&self) -> RT {
"
        Return lists of orbital in mixed spin index alternating alpha and beta
        ";
return Determinant::obtBits2ObtMixSpinIndexList(self.alphaObtBits, self.betaObtBits);
}
fn countNumOrbitalsInBits<T0, RT>(bits: T0) -> RT {
"
        Return the number of orbitals in this bits
        ";
let mut count = 0;
while bits != 0 {
if (bits & 1) == 1 {
count += 1;
}
bits >>= 1;
}
return count;
}
fn countNumOrbitalsInBitsUpTo4<T0, RT>(bits: T0) -> RT {
"
        Return the number of orbitals in this bits
        ";
let mut count = 0;
while bits != 0&&count < 4 {
if (bits & 1) == 1 {
count += 1;
}
bits >>= 1;
}
return count;
}
fn obtBits2ObtIndexList<T0, RT>(bits: T0) -> RT {
"
        Return the corresponding list of orbital numbers from orbital bits
        ";
let mut i = 0;
let mut obts = vec![];
while bits != 0 {
if (bits & 1) == 1 {
obts.push(i);
}
bits >>= 1;
i += 1;
}
return obts;
}
fn mixIndexList<T0, T1, RT>(alphaList: T0, betaList: T1) -> RT {
"
        Mix the alpha and beta orbital index list to one mixed list
        ";
return (alphaList.iter().map(|elem| (elem*2)).collect::<Vec<_>>() + betaList.iter().map(|elem| ((elem*2) + 1)).collect::<Vec<_>>());
}
fn obtBits2ObtMixSpinIndexList<T0, T1, RT>(alphaBits: T0, betaBits: T1) -> RT {
"
        Return the corresponding list of orbital numbers of orbital bits
        ";
let (alphaList, betaList) = (Determinant::obtBits2ObtIndexList(alphaBits), Determinant::obtBits2ObtIndexList(betaBits));
return Determinant::mixIndexList(alphaList, betaList);
}
fn obtIndexList2ObtBits<T0, RT>(obtList: T0) -> RT {
"
        Return the corresponding orbital bits of list from orbital numbers
        ";
if obtList.len() == 0 {
return 0;
}
obtList = sorted(obtList, true);
let mut iPre = obtList[0];
let mut bits = 1;
for i in obtList {
bits <<= (iPre - i);
bits |= 1;
iPre = i;
}
bits <<= iPre;
return bits;
}
fn getOrbitalPositions<T0, T1, RT>(bits: T0, orbitalIndexList: T1) -> RT {
"
        Return the position of orbital in determinant
        ";
let mut count = 0;
let mut index = 0;
let mut positions = vec![];
for i in orbitalIndexList {
while index < i {
if (bits & 1) == 1 {
count += 1;
}
bits >>= 1;
index += 1;
}
positions.push(count);
continue;
}
return positions;
}
fn getOrbitalPositionLists<T0, T1, RT>(&self, alphaIndexList: T0, betaIndexList: T1) -> RT {
"
        Return the positions of indexes in lists
        ";
return (Determinant::getOrbitalPositions(self.alphaObtBits, alphaIndexList), Determinant::getOrbitalPositions(self.betaObtBits, betaIndexList));
}
fn addAlphaOrbital<T0>(&self, orbitalIndex: T0)  {
"
        Add an alpha orbital to the determinant
        ";
self.alphaObtBits |= (1 << orbitalIndex);
}
fn addBetaOrbital<T0>(&self, orbitalIndex: T0)  {
"
        Add an beta orbital to the determinant
        ";
self.betaObtBits |= (1 << orbitalIndex);
}
fn removeAlphaOrbital<T0>(&self, orbitalIndex: T0)  {
"
        Remove an alpha orbital from the determinant
        ";
self.alphaObtBits &= None(1 << orbitalIndex);
}
fn removeBetaOrbital<T0>(&self, orbitalIndex: T0)  {
"
        Remove an beta orbital from the determinant
        ";
self.betaObtBits &= None(1 << orbitalIndex);
}
fn numberOfCommonOrbitals<T0, RT>(&self, another: T0) -> RT {
"
        Return the number of common orbitals between this determinant and another determinant
        ";
return (Determinant::countNumOrbitalsInBits((self.alphaObtBits & another.alphaObtBits)), Determinant::countNumOrbitalsInBits((self.betaObtBits & another.betaObtBits)));
}
fn getCommonOrbitalsInLists<T0, RT>(&self, another: T0) -> RT {
"Return common orbitals between this determinant and another determinant in lists";
return (Determinant::obtBits2ObtIndexList((self.alphaObtBits & another.alphaObtBits)), Determinant::obtBits2ObtIndexList((self.betaObtBits & another.betaObtBits)));
}
fn getCommonOrbitalsInMixedSpinIndexList<T0, RT>(&self, another: T0) -> RT {
let (alphaList, betaList) = self.getCommonOrbitalsInLists(another);
return Determinant::mixIndexList(alphaList, betaList);
}
fn numberOfDiffOrbitals<T0, RT>(&self, another: T0) -> RT {
"
        Return the number of different alpha and beta orbitals between this determinant and another determinant
        ";
let (diffAlpha, diffBeta) = (Determinant::countNumOrbitalsInBits((self.alphaObtBits ^ another.alphaObtBits)), Determinant::countNumOrbitalsInBits((self.betaObtBits ^ another.betaObtBits)));
return ((diffAlpha/2), (diffBeta/2));
}
fn numberOfTotalDiffOrbitals<T0, RT>(&self, another: T0) -> RT {
"
        Return the number of different orbitals between this determinant and another determinant
        ";
let (diffAlpha, diffBeta) = self.numberOfDiffOrbitals(another);
return (diffAlpha + diffBeta);
}
fn diff2OrLessOrbitals<T0, RT>(&self, another: T0) -> RT {
"
        Return true if two determinants differ 2 or less orbitals
        ";
let (diffAlpha, diffBeta) = (Determinant::countNumOrbitalsInBitsUpTo4((self.alphaObtBits ^ another.alphaObtBits)), Determinant::countNumOrbitalsInBitsUpTo4((self.betaObtBits ^ another.betaObtBits)));
return (diffAlpha + diffBeta) <= 4;
}
fn uniqueOrbitalsInBits<T0, T1, RT>(bits1: T0, bits2: T1) -> RT {
"
        Return the unique bits in two different bits
        ";
let common = (bits1 & bits2);
return ((bits1 ^ common), (bits2 ^ common));
}
fn uniqueOrbitalsInLists<T0, T1, RT>(bits1: T0, bits2: T1) -> RT {
"
        Return the unique bits in two different bits
        ";
let (bits1, bits2) = Determinant::uniqueOrbitalsInBits(bits1, bits2);
return (Determinant::obtBits2ObtIndexList(bits1), Determinant::obtBits2ObtIndexList(bits2));
}
fn getUniqueOrbitalsInLists<T0, RT>(&self, another: T0) -> RT {
"
        Return the unique orbital lists in two different determinants
        ";
let (alphaList1, alphaList2) = Determinant::uniqueOrbitalsInLists(self.alphaObtBits, another.alphaObtBits);
let (betaList1, betaList2) = Determinant::uniqueOrbitalsInLists(self.betaObtBits, another.betaObtBits);
return ((alphaList1, betaList1), (alphaList2, betaList2));
}
fn getUnoccupiedOrbitalsInLists<T0, RT>(&self, nmo: T0) -> RT {
"
        Return the unoccupied orbitals in the determinants
        ";
let mut alphaBits = Noneself.alphaObtBits;
let mut betaBits = Noneself.betaObtBits;
let mut alphaObts = vec![];
let mut betaObts = vec![];
for i in (0..nmo) {
if (alphaBits & 1) == 1 {
alphaObts.push(i);
}
alphaBits >>= 1;
if (betaBits & 1) == 1 {
betaObts.push(i);
}
betaBits >>= 1;
}
return (alphaObts, betaObts);
}
fn getSignToMoveOrbitalsToTheFront<T0, T1, RT>(&self, alphaIndexList: T0, betaIndexList: T1) -> RT {
"
        Return the final sign if move listed orbitals to the front
        ";
let mut sign = 1;
let (alphaPositions, betaPositions) = self.getOrbitalPositionLists(alphaIndexList, betaIndexList);
for i in (0..alphaPositions.len()) {
if ((alphaPositions[i] - i) % 2) == 1 {
sign = -(sign);
}
}
for i in (0..betaPositions.len()) {
if ((betaPositions[i] - i) % 2) == 1 {
sign = -(sign);
}
}
return sign;
}
fn getUniqueOrbitalsInListsPlusSign<T0, RT>(&self, another: T0) -> RT {
"
        Return the unique orbital lists in two different determinants and the sign of the maximum coincidence determinants
        ";
let (alphaList1, alphaList2) = Determinant::uniqueOrbitalsInLists(self.alphaObtBits, another.alphaObtBits);
let (betaList1, betaList2) = Determinant::uniqueOrbitalsInLists(self.betaObtBits, another.betaObtBits);
let (sign1, sign2) = (self.getSignToMoveOrbitalsToTheFront(alphaList1, betaList1), another.getSignToMoveOrbitalsToTheFront(alphaList2, betaList2));
return ((alphaList1, betaList1), (alphaList2, betaList2), (sign1*sign2));
}
fn getUniqueOrbitalsInMixIndexListsPlusSign<T0, RT>(&self, another: T0) -> RT {
"
        Return the unique orbital lists in two different determinants and the sign of the maximum coincidence determinants
        ";
let ((alphaList1, betaList1), (alphaList2, betaList2), sign) = self.getUniqueOrbitalsInListsPlusSign(another);
return (Determinant::mixIndexList(alphaList1, betaList1), Determinant::mixIndexList(alphaList2, betaList2), sign);
}
fn toIntTuple<RT>(&self) -> RT {
"
        Return a int tuple
        ";
return (self.alphaObtBits, self.betaObtBits);
}
fn createFromIntTuple<T0, RT>(intTuple: T0) -> RT {
return Determinant(intTuple[0], intTuple[1]);
}
fn generateSingleExcitationsOfDet<T0, RT>(&self, nmo: T0) -> RT {
"
        Generate all the single excitations of determinant in a list
        ";
let (alphaO, betaO) = self.getOrbitalIndexLists();
let (alphaU, betaU) = self.getUnoccupiedOrbitalsInLists(nmo);
let mut dets = vec![];
for i in alphaO {
for j in alphaU {
let mut det = self.copy();
det.removeAlphaOrbital(i);
det.addAlphaOrbital(j);
dets.push(det);
}
}
for k in betaO {
for l in betaU {
let mut det = self.copy();
det.removeBetaOrbital(k);
det.addBetaOrbital(l);
dets.push(det);
}
}
return dets;
}
fn generateDoubleExcitationsOfDet<T0, RT>(&self, nmo: T0) -> RT {
"
        Generate all the double excitations of determinant in a list
        ";
let (alphaO, betaO) = self.getOrbitalIndexLists();
let (alphaU, betaU) = self.getUnoccupiedOrbitalsInLists(nmo);
let mut dets = vec![];
for i in alphaO {
for j in alphaU {
for k in betaO {
for l in betaU {
let mut det = self.copy();
det.removeAlphaOrbital(i);
det.addAlphaOrbital(j);
det.removeBetaOrbital(k);
det.addBetaOrbital(l);
dets.push(det);
}
}
}
}
for (i1, i2) in combinations(alphaO, 2) {
for (j1, j2) in combinations(alphaU, 2) {
let mut det = self.copy();
det.removeAlphaOrbital(i1);
det.addAlphaOrbital(j1);
det.removeAlphaOrbital(i2);
det.addAlphaOrbital(j2);
dets.push(det);
}
}
for (k1, k2) in combinations(betaO, 2) {
for (l1, l2) in combinations(betaU, 2) {
let mut det = self.copy();
det.removeBetaOrbital(k1);
det.addBetaOrbital(l1);
det.removeBetaOrbital(k2);
det.addBetaOrbital(l2);
dets.push(det);
}
}
return dets;
}
fn generateSingleAndDoubleExcitationsOfDet<T0, RT>(&self, nmo: T0) -> RT {
"
        Generate all the single and double excitations of determinant in a list
        ";
return (self.generateSingleExcitationsOfDet(nmo) + self.generateDoubleExcitationsOfDet(nmo));
}
fn copy<RT>(&self) -> RT {
"
        Return a deep copy of self
        ";
return Determinant(self.alphaObtBits, self.betaObtBits);
}
fn __str__<RT>(&self) -> RT {
"
        Print a representation of the Determinant
        ";
let (a, b) = self.getOrbitalIndexLists();
return ((("|" + String::from(a)) + String::from(b)) + ">");
} 
}
struct HamiltonianGenerator {
Hspin: ST0,
antiSym2eInt: ST1,
}

impl HamiltonianGenerator {
"
    class for Full CI matrix elements
    ";
fn __init__<T0, T1>(&self, H_spin: T0, mo_spin_eri: T1)  {
"
        Constructor for MatrixElements
        ";
self.Hspin = H_spin;
self.antiSym2eInt = mo_spin_eri;
}
fn generateMatrix<T0, RT>(&self, detList: T0) -> RT {
"
        Generate CI Matrix
        ";
let numDet = detList.len();
let matrix = np.zeros((numDet, numDet));
for i in (0..numDet) {
for j in (0..(i + 1)) {
matrix[(i, j)] = self.calcMatrixElement(detList[i], detList[j]);
matrix[(j, i)] = matrix[(i, j)];
}
}
return matrix;
}
fn calcMatrixElement<T0, T1, RT>(&self, det1: T0, det2: T1) -> RT {
"
        Calculate a matrix element by two determinants
        ";
let mut numUniqueOrbitals = None;
if det1.diff2OrLessOrbitals(det2) {
numUniqueOrbitals = det1.numberOfTotalDiffOrbitals(det2);
if numUniqueOrbitals == 0 {
return self.calcMatrixElementIdentialDet(det1);
}
if numUniqueOrbitals == 2 {
return self.calcMatrixElementDiffIn2(det1, det2);
} else {
if numUniqueOrbitals == 1 {
return self.calcMatrixElementDiffIn1(det1, det2);
} else {
return 0.0;
}
}
} else {
return 0.0;
}
}
fn calcMatrixElementDiffIn2<T0, T1, RT>(&self, det1: T0, det2: T1) -> RT {
"
        Calculate a matrix element by two determinants where the determinants differ by 2 spin orbitals
        ";
let (unique1, unique2, sign) = det1.getUniqueOrbitalsInMixIndexListsPlusSign(det2);
return (sign*self.antiSym2eInt[(unique1[0], unique1[1], unique2[0], unique2[1])]);
}
fn calcMatrixElementDiffIn1<T0, T1, RT>(&self, det1: T0, det2: T1) -> RT {
"
        Calculate a matrix element by two determinants where the determinants differ by 1 spin orbitals
        ";
let (unique1, unique2, sign) = det1.getUniqueOrbitalsInMixIndexListsPlusSign(det2);
let m = unique1[0];
let p = unique2[0];
let Helem = self.Hspin[(m, p)];
let common = det1.getCommonOrbitalsInMixedSpinIndexList(det2);
let mut Relem = 0.0;
for n in common {
Relem += self.antiSym2eInt[(m, n, p, n)];
}
return (sign*(Helem + Relem));
}
fn calcMatrixElementIdentialDet<T0, RT>(&self, det: T0) -> RT {
"
        Calculate a matrix element by two determinants where they are identical
        ";
let spinObtList = det.getOrbitalMixedIndexList();
let mut Helem = 0.0;
for m in spinObtList {
Helem += self.Hspin[(m, m)];
}
let length = spinObtList.len();
let mut Relem = 0.0;
for m in (0..(length - 1)) {
for n in ((m + 1)..length) {
Relem += self.antiSym2eInt[(spinObtList[m], spinObtList[n], spinObtList[m], spinObtList[n])];
}
}
return (Helem + Relem);
} 
}