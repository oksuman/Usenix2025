package ApproxSS

import (
	"fmt"
	"io"
	"math/big"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v4/bfv"
	"github.com/tuneinsight/lattigo/v4/drlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type MyShamirApproxSS struct {
	VanSS *VanillaShamirSS

	params4DoubleEncryption bfv.Parameters
	params4cmb              rlwe.Parameters
	f                       *skEncryptionBFV
	ek1All                  []*rlwe.SecretKey
	ek2All                  []*rlwe.SecretKey
	N                       int
	T                       int

	// Added: store last measured times for each stage
	LastR1Party time.Duration
	LastR1Agg   time.Duration
	LastR2Party time.Duration
	LastR2Agg   time.Duration
}

func NewMyShamirApproxSS(N, T int, params bfv.Parameters) (myShamirApproxSS *MyShamirApproxSS) {
	myShamirApproxSS = new(MyShamirApproxSS)
	myShamirApproxSS.N = N
	myShamirApproxSS.T = T
	myShamirApproxSS.VanSS = NewVanillaShamirSS(N, T, params)
	myShamirApproxSS.params4DoubleEncryption = findLargerParameters(params, N)

	myShamirApproxSS.f = NewSKencryption(myShamirApproxSS.params4DoubleEncryption)
	myShamirApproxSS.params4cmb, _ = rlwe.NewParameters(params.LogN(), params.Q()[:1], nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), params.DefaultNTTFlag())

	myShamirApproxSS.ek1All = make([]*rlwe.SecretKey, N)
	myShamirApproxSS.ek2All = make([]*rlwe.SecretKey, N)
	return
}

func findLargerParameters(params bfv.Parameters, N int) (params4DoubleEncryption bfv.Parameters) {
	testTimesforFEparams := 10000
	indexForQ1 := 1

	T := params.Q()[0]
	Tbig := new(big.Int).SetUint64(T)
	bound := new(big.Int).Mul(big.NewInt(int64(2*N*int(params.NoiseBound()))), Tbig)

	newQ := make([]uint64, 2)
	newQ[0] = T
	for i := int(T*uint64(N)*params.NoiseBound()) / params.N(); i < int(T*uint64(N)*params.NoiseBound())/params.N()+testTimesforFEparams; i++ {
		x := int(T) + i*2*params.N()
		if params.RingType() == ring.ConjugateInvariant {
			x = x + i*2*params.N()
		}
		xBig := big.NewInt(int64(x))
		if xBig.Cmp(bound) == 1 {
			if xBig.ProbablyPrime(0) {
				indexForQ1 = indexForQ1 - 1
				if indexForQ1 == 0 {
					newQ[1] = uint64(x)
					FEparams, err := rlwe.NewParameters(params.LogN(), newQ, nil, 0, params.HammingWeight(), params.Sigma(), params.RingType(), rlwe.NewScale(1), false)
					if err != nil {
						panic(err)
					}
					bfvParams, _ := bfv.NewParameters(FEparams, FEparams.Q()[0])
					return bfvParams
				}
			}
		}
	}
	panic("The number of loops is not large enough to find proper parameters of double encryption!")
}

func (myShamirApproxSS *MyShamirApproxSS) ApproxRecover4TestTime(bound4smudgingNoise int) (isSuc bool, timeCompOnce time.Duration, timeCompNotOnce time.Duration, sizeCommOnce float64, sizeCommNotOnce float64) {

	paramsMesSpace := myShamirApproxSS.VanSS.thdizer.params
	ringQ := paramsMesSpace.RingQ()
	levelQ := 0

	var timeStage1, timeStage2, timeStage3, timeStage4 time.Duration
	var timeNotOnce1, timeNotOnce2 time.Duration

	// -------------------------
	// Round 1 starts
	// -------------------------
	CRS, _ := utils.NewPRNG()
	a := TagGen(myShamirApproxSS.params4DoubleEncryption, CRS)
	parR1 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)

	CTni_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	CTsi_all := make(map[int]*rlwe.Ciphertext, myShamirApproxSS.VanSS.T)
	lagrangeCoeffs := make(map[int]uint64, myShamirApproxSS.VanSS.T)

	// Round1 - Party (simulate: compute once, reuse)
	for i, par := range parR1 {
		if i == 0 {
			filename := fmt.Sprintf("%s%d", "mesToParty", par)
			skFilename := fmt.Sprintf("%s%s", "file/skShare/", filename)
			skFile, _ := os.Open(skFilename)
			skData, _ := io.ReadAll(skFile)
			skFile.Close()

			share := myShamirApproxSS.VanSS.thdizer.AllocateThresholdSecretShare()
			share.UnmarshalBinary(skData)

			myShamirApproxSS.ek1All[par], timeStage1 = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "s")
			myShamirApproxSS.ek2All[par], _ = myShamirApproxSS.f.GenerateEKThenShareToFile(myShamirApproxSS.N, myShamirApproxSS.T, par, "n")
			timeStage1 = timeStage1 * 2

			timeNotOnce1Start := time.Now()
			prng_i, _ := utils.NewPRNG()
			smudgingNoiseSampler := ring.NewGaussianSampler(prng_i, paramsMesSpace.RingQ(), float64(bound4smudgingNoise)/6, bound4smudgingNoise)
			ni := ringQ.NewPoly()
			smudgingNoiseSampler.ReadLvl(levelQ, ni)

			CTsi_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek1All[par], share.Poly.Q)
			CTni_all[par] = myShamirApproxSS.f.EncPolyCoeff(a, myShamirApproxSS.ek2All[par], ni)
			timeNotOnce1 = time.Since(timeNotOnce1Start)
		} else {
			myShamirApproxSS.ek1All[par] = myShamirApproxSS.ek1All[parR1[0]]
			myShamirApproxSS.ek2All[par] = myShamirApproxSS.ek2All[parR1[0]]
			CTsi_all[par] = CTsi_all[parR1[0]]
			CTni_all[par] = CTni_all[parR1[0]]
		}
	}

	// Round1 - Aggregator
	timeStage2Start := time.Now()
	survivialPublicPoint := []drlwe.ShamirPublicPoint{}
	for _, pi := range parR1 {
		survivialPublicPoint = append(survivialPublicPoint, drlwe.ShamirPublicPoint(pi+1))
	}
	cmbR1 := NewMyCombiner(&myShamirApproxSS.params4cmb, survivialPublicPoint, myShamirApproxSS.VanSS.T)
	for _, cIdx := range parR1 {
		lagrangeCoeffs[cIdx] = fromRNSScalarToInt(myShamirApproxSS.params4cmb.RingQP(),
			cmbR1.LagrangeCoeffs[drlwe.ShamirPublicPoint(cIdx+1)])
	}
	timeStage2 = time.Since(timeStage2Start)

	// -------------------------
	// Round 2 starts
	// -------------------------
	parR2 := SelectRandomParticipants(myShamirApproxSS.VanSS.N, myShamirApproxSS.VanSS.T)
	dkShare_All := make(map[int]drlwe.ShamirSecretShare, myShamirApproxSS.VanSS.T)

	// Round2 - Party
	for i, par := range parR2 {
		if i == 0 {
			dkShare_All[par], timeStage3 = myShamirApproxSS.f.GenerateDKShareFile4TestTime(parR1, lagrangeCoeffs, par)
		} else {
			dkShare_All[par] = dkShare_All[parR2[0]]
		}
	}

	// Round2 - Aggregator
	timeStage4Start := time.Now()
	DK := myShamirApproxSS.f.GenerateDK(parR2, dkShare_All)
	timeStage4 = time.Since(timeStage4Start)

	timeNotOnce2Start := time.Now()
	approxMessage := ringQ.NewPoly()
	myShamirApproxSS.f.FEDecFinal(DK, CTni_all, CTsi_all, lagrangeCoeffs, parR1, approxMessage)
	timeNotOnce2 = time.Since(timeNotOnce2Start)

	// -------------------------
	// Summary and reporting
	// -------------------------
	timeRound1Party := timeStage1
	timeRound1Agg := timeStage2
	timeRound2Party := timeStage3
	timeRound2Agg := timeStage4
	timeCompTotal := timeRound1Party + timeRound1Agg + timeRound2Party + timeRound2Agg

	// Store for main.go access
	myShamirApproxSS.LastR1Party = timeRound1Party
	myShamirApproxSS.LastR1Agg = timeRound1Agg
	myShamirApproxSS.LastR2Party = timeRound2Party
	myShamirApproxSS.LastR2Agg = timeRound2Agg

	// Console output
	fmt.Println("========================================")
	fmt.Printf("N=%d  T=%d  (Computation-only reporting)\n", myShamirApproxSS.N, myShamirApproxSS.T)
	fmt.Println("----------------------------------------")
	fmt.Printf("Round1-Party : %v\n", timeRound1Party)
	fmt.Printf("Round1-Agg   : %v\n", timeRound1Agg)
	fmt.Printf("Round2-Party : %v\n", timeRound2Party)
	fmt.Printf("Round2-Agg   : %v\n", timeRound2Agg)
	fmt.Println("----------------------------------------")
	fmt.Printf("Total-Comp   : %v\n", timeCompTotal)
	fmt.Println("========================================")

	timeCompOnce = timeCompTotal
	timeCompNotOnce = timeNotOnce1 + timeNotOnce2
	sizeCommOnce = 0.0
	sizeCommNotOnce = 0.0

	isSuc = TestResult(approxMessage, myShamirApproxSS.VanSS.secret.Value.Q,
		uint64(myShamirApproxSS.VanSS.T*bound4smudgingNoise),
		myShamirApproxSS.VanSS.thdizer.params.Q()[0:1])
	return
}
