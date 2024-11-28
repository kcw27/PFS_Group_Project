package main

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"time"
)

// Helper function to calculate Spearman correlation
func spearmanCorrelation(x, y []float64) float64 {
	if len(x) != len(y) {
		panic("Input slices must have equal length")
	}

	n := len(x)
	
	// Create rank pairs
	type pair struct {
		value float64
		rank  float64
	}
	
	// Calculate ranks for x
	xPairs := make([]pair, n)
	for i, v := range x {
		xPairs[i] = pair{v, float64(i)}
	}
	sort.Slice(xPairs, func(i, j int) bool {
		return xPairs[i].value < xPairs[j].value
	})
	for i := range xPairs {
		xPairs[i].rank = float64(i + 1)
	}
	
	// Calculate ranks for y
	yPairs := make([]pair, n)
	for i, v := range y {
		yPairs[i] = pair{v, float64(i)}
	}
	sort.Slice(yPairs, func(i, j int) bool {
		return yPairs[i].value < yPairs[j].value
	})
	for i := range yPairs {
		yPairs[i].rank = float64(i + 1)
	}
	
	// Restore original order and calculate correlation
	xRanks := make([]float64, n)
	yRanks := make([]float64, n)
	for i := range xPairs {
		idx := int(xPairs[i].rank) - 1
		xRanks[idx] = float64(i + 1)
	}
	for i := range yPairs {
		idx := int(yPairs[i].rank) - 1
		yRanks[idx] = float64(i + 1)
	}
	
	// Calculate mean of ranks
	var sumX, sumY float64
	for i := 0; i < n; i++ {
		sumX += xRanks[i]
		sumY += yRanks[i]
	}
	meanX := sumX / float64(n)
	meanY := sumY / float64(n)
	
	// Calculate correlation
	var numerator, denomX, denomY float64
	for i := 0; i < n; i++ {
		dx := xRanks[i] - meanX
		dy := yRanks[i] - meanY
		numerator += dx * dy
		denomX += dx * dx
		denomY += dy * dy
	}
	
	return numerator / math.Sqrt(denomX*denomY)
}

// Helper function to get columns matching a specific color
func getColumnsForColor(data [][]float64, colorh1C1C2 map[string]string, color string) [][]float64 {
	var columns [][]float64
	for i := range data[0] {
		if val, exists := colorh1C1C2[fmt.Sprintf("gene%d", i+1)]; exists && val == color {
			column := make([]float64, len(data))
			for j := range data {
				column[j] = data[j][i]
			}
			columns = append(columns, column)
		}
	}
	return columns
}

// dispersionModule2Module calculates the dispersion value between two modules
func dispersionModule2Module(c1, c2 string, datC1, datC2 [][]float64, colorh1C1C2 map[string]string) float64 {
	if c1 == c2 {
		// Get columns for the module
		columnsC1 := getColumnsForColor(datC1, colorh1C1C2, c1)
		columnsC2 := getColumnsForColor(datC2, colorh1C1C2, c1)
		
		n := len(columnsC1)
		if n == 0 {
			return 0.0
		}
		
		// Calculate correlation matrices
		var sumDifCorSquared float64
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				corC1 := spearmanCorrelation(columnsC1[i], columnsC1[j])
				corC2 := spearmanCorrelation(columnsC2[i], columnsC2[j])
				difCor := corC1 - corC2
				sumDifCorSquared += difCor * difCor
			}
		}
		
		// Calculate final dispersion
		denominator := float64(n*n-n) / 2.0
		return math.Sqrt((1.0 / denominator) * (sumDifCorSquared / 2.0))
		
	} else {
		// Get columns for both modules
		columnsC1_1 := getColumnsForColor(datC1, colorh1C1C2, c1)
		columnsC1_2 := getColumnsForColor(datC1, colorh1C1C2, c2)
		columnsC2_1 := getColumnsForColor(datC2, colorh1C1C2, c1)
		columnsC2_2 := getColumnsForColor(datC2, colorh1C1C2, c2)
		
		n1 := len(columnsC1_1)
		n2 := len(columnsC1_2)
		if n1 == 0 || n2 == 0 {
			return 0.0
		}
		
		// Calculate correlation differences
		var sumDifCorSquared float64
		for i := 0; i < n1; i++ {
			for j := 0; j < n2; j++ {
				corC1 := spearmanCorrelation(columnsC1_1[i], columnsC1_2[j])
				corC2 := spearmanCorrelation(columnsC2_1[i], columnsC2_2[j])
				difCor := corC1 - corC2
				sumDifCorSquared += difCor * difCor
			}
		}
		
		// Calculate final dispersion
		return math.Sqrt((1.0 / float64(n1*n2)) * sumDifCorSquared)
	}
}

// generatePermutations creates a set of permuted indices
func generatePermutations(datC1, datC2 [][]float64, numPermutations int) [][]int {
	totalSamples := len(datC1) + len(datC2)
	sampleSize := len(datC1)
	
	// Initialize the permutations slice
	permutations := make([][]int, numPermutations)
	
	// Create random number generator with time-based seed
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	
	// Generate permutations
	for i := 0; i < numPermutations; i++ {
		// Create a slice of all possible indices
		indices := make([]int, totalSamples)
		for j := range indices {
			indices[j] = j
		}
		
		// Shuffle the indices using Fisher-Yates algorithm
		for j := totalSamples - 1; j > 0; j-- {
			k := r.Intn(j + 1)
			indices[j], indices[k] = indices[k], indices[j]
		}
		
		// Take the first sampleSize elements as our permutation
		permutations[i] = make([]int, sampleSize)
		copy(permutations[i], indices[:sampleSize])
	}
	
	return permutations
}

// scaleData scales the input data to have mean 0 and variance 1
func scaleData(data [][]float64) [][]float64 {
    rows := len(data)
    if rows == 0 {
        return data
    }
    cols := len(data[0])
    
    // Create scaled data matrix
    scaledData := make([][]float64, rows)
    for i := range scaledData {
        scaledData[i] = make([]float64, cols)
    }
    
    // Scale each column
    for j := 0; j < cols; j++ {
        // Calculate mean
        var sum float64
        for i := 0; i < rows; i++ {
            sum += data[i][j]
        }
        mean := sum / float64(rows)
        
        // Calculate variance
        var sumSquaredDiff float64
        for i := 0; i < rows; i++ {
            diff := data[i][j] - mean
            sumSquaredDiff += diff * diff
        }
        stdDev := math.Sqrt(sumSquaredDiff / float64(rows-1))
        
        // Scale the column
        if stdDev != 0 {
            for i := 0; i < rows; i++ {
                scaledData[i][j] = (data[i][j] - mean) / stdDev
            }
        } else {
            // If stdDev is 0, just center the data
            for i := 0; i < rows; i++ {
                scaledData[i][j] = data[i][j] - mean
            }
        }
    }
    
    return scaledData
}

// combineAndScaleData combines and scales datC1 and datC2
func combineAndScaleData(datC1, datC2 [][]float64) [][]float64 {
    // Calculate total rows needed
    totalRows := len(datC1) + len(datC2)
    if totalRows == 0 {
        return nil
    }
    
    // Create combined data matrix
    combinedData := make([][]float64, totalRows)
    
    // Copy datC1
    for i := range datC1 {
        combinedData[i] = make([]float64, len(datC1[i]))
        copy(combinedData[i], datC1[i])
    }
    
    // Copy datC2
    for i := range datC2 {
        combinedData[i+len(datC1)] = make([]float64, len(datC2[i]))
        copy(combinedData[i+len(datC1)], datC2[i])
    }
    
    // Scale the combined data
    return scaleData(combinedData)
}

// permutationProcedureModule2Module calculates dispersion values using permuted data
func permutationProcedureModule2Module(permutation []int, d [][]float64, c1, c2 string, colorh1C1C2 map[string]string) float64 {
    // Create d1 from permuted indices
    d1 := make([][]float64, len(permutation))
    for i, idx := range permutation {
        d1[i] = make([]float64, len(d[0]))
        copy(d1[i], d[idx])
    }
    
    // Create d2 from remaining indices (complement of permutation)
    d2 := make([][]float64, len(d)-len(permutation))
    usedIndices := make(map[int]bool)
    for _, idx := range permutation {
        usedIndices[idx] = true
    }
    
    j := 0
    for i := 0; i < len(d); i++ {
        if !usedIndices[i] {
            d2[j] = make([]float64, len(d[0]))
            copy(d2[j], d[i])
            j++
        }
    }
    
    // Calculate dispersion using permuted datasets
    return dispersionModule2Module(c1, c2, d1, d2, colorh1C1C2)
}


/*

package main

import (
	"fmt"
)

func main() {
	// Example data
	datC1 := [][]float64{{1, 2}, {3, 4}}
	datC2 := [][]float64{{5, 6}, {7, 8}}
	colorh1C1C2 := map[string]string{"gene1": "red", "gene2": "blue"}

	// Generate permutations
	numPermutations := 1000
	permutations := generatePermutations(datC1, datC2, numPermutations)

	// Scale and combine the data
	d := combineAndScaleData(datC1, datC2)

	// Compute dispersion matrix and null distribution
	uniqueColors := []string{"red", "blue"} // Example unique colors
	dispersionMatrix := make([][]float64, len(uniqueColors))
	nullDistrib := make(map[string]map[string][]float64)

	for i, c1 := range uniqueColors {
		dispersionMatrix[i] = make([]float64, len(uniqueColors))
		nullDistrib[c1] = make(map[string][]float64)
		for j, c2 := range uniqueColors {
			// Calculate dispersion
			dispersionMatrix[i][j] = dispersionModule2Module(c1, c2, datC1, datC2, colorh1C1C2)

			// Generate null distribution
			nullDistrib[c1][c2] = make([]float64, numPermutations)
			for k, perm := range permutations {
				nullDistrib[c1][c2][k] = permutationProcedureModule2Module(perm, d, c1, c2, colorh1C1C2)
			}
		}
	}

	// Create permutation summary
	permutationSummary := make([][]int, len(uniqueColors))
	for i := range permutationSummary {
		permutationSummary[i] = make([]int, len(uniqueColors))
		for j := range permutationSummary[i] {
			count := 0
			for _, val := range nullDistrib[uniqueColors[i]][uniqueColors[j]] {
				if val >= dispersionMatrix[i][j] {
					count++
				}
			}
			permutationSummary[i][j] = count
		}
	}

	// Print results
	fmt.Println("Dispersion Matrix:", dispersionMatrix)
	fmt.Println("Permutation Summary:", permutationSummary)
}

*/