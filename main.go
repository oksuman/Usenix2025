package main

import (
	"ApproxSS/ApproxSS"
	"bufio"
	"flag"
	"fmt"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v4/bfv"
)

type ParamPair struct {
	N int
	T int
}

func main() {
	var B = flag.Int("B", 16, "The bound of each nonce")
	flag.Parse()

	params, err := bfv.NewParametersFromLiteral(bfv.PN12QP109)
	if err != nil {
		fmt.Println(err)
		return
	}

	paramPairs := []ParamPair{
		{32, 5}, {32, 10}, {32, 16}, {32, 21},
		{64, 5}, {64, 6}, {64, 21}, {64, 32}, {64, 42},
		{128, 5}, {128, 7}, {128, 42}, {128, 64}, {128, 85},
		{256, 5}, {256, 8}, {256, 85}, {256, 128}, {256, 170},
		{512, 5}, {512, 9}, {512, 170}, {512, 256}, {512, 341},
		{1024, 5}, {1024, 10}, {1024, 341}, {1024, 512}, {1024, 682},
	}

	const trials = 1
	filename := "TimeComp.txt"

	for _, pair := range paramPairs {
		n := pair.N
		t := pair.T

		var sumR1Party, sumR1Agg, sumR2Party, sumR2Agg, sumTotal time.Duration

		for i := 0; i < trials; i++ {
			myShamir := ApproxSS.NewMyShamirApproxSS(n, t, params)
			myShamir.VanSS.ShareThenWrite(nil, "skShare", nil)

			_, timeCompOnce, _, _, _ := myShamir.ApproxRecover4TestTime(*B)

			// Extract the per-stage times from ApproxRecover4TestTime console output if needed
			sumR1Party += myShamir.LastR1Party
			sumR1Agg += myShamir.LastR1Agg
			sumR2Party += myShamir.LastR2Party
			sumR2Agg += myShamir.LastR2Agg
			sumTotal += timeCompOnce
		}

		avgR1Party := sumR1Party / trials
		avgR1Agg := sumR1Agg / trials
		avgR2Party := sumR2Party / trials
		avgR2Agg := sumR2Agg / trials
		avgTotal := sumTotal / trials

		ringDim := params.N()
		logQ := params.LogQP()

		RecordAverageDetailed(filename, n, t, ringDim, logQ,
			avgR1Party, avgR1Agg, avgR2Party, avgR2Agg, avgTotal)

		fmt.Printf("N=%d, T=%d | AvgTotal: %.6fs\n", n, t, avgTotal.Seconds())
	}
}

// RecordAverageDetailed writes averaged per-round times + BFV params
func RecordAverageDetailed(filename string, N, T, ringDim, logQ int,
	avgR1Party, avgR1Agg, avgR2Party, avgR2Agg, avgTotal time.Duration) {

	f, err := os.OpenFile(filename, os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
	if err != nil {
		fmt.Println("Error opening file:", err)
		return
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	now := time.Now().Format("2006-01-02 15:04:05")

	w.WriteString(fmt.Sprintf("=============================================\n"))
	w.WriteString(fmt.Sprintf("ATASSES Performance Report (%s)\n", now))
	w.WriteString(fmt.Sprintf("Parameters: N = %d, T = %d\n", N, T))
	w.WriteString(fmt.Sprintf("BFV Params: RingDim = %d, logQ = %d bits\n", ringDim, logQ))
	w.WriteString("---------------------------------------------\n")
	w.WriteString(fmt.Sprintf("Avg Round1-Party : %v\n", avgR1Party))
	w.WriteString(fmt.Sprintf("Avg Round1-Agg   : %v\n", avgR1Agg))
	w.WriteString(fmt.Sprintf("Avg Round2-Party : %v\n", avgR2Party))
	w.WriteString(fmt.Sprintf("Avg Round2-Agg   : %v\n", avgR2Agg))
	w.WriteString(fmt.Sprintf("Avg Total        : %v\n", avgTotal))
	w.WriteString(fmt.Sprintf("Avg Total (sec)  : %.6f\n", avgTotal.Seconds()))
	w.WriteString("=============================================\n\n")
	w.Flush()
}
