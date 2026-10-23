//! Selection of the out-going cell at a fracture intersection: complete
//! mixing and streamline routing (DFNTrans `CompleteMixingRandomSampling`,
//! `StreamlineRandomSampling`).

use crate::random::Drand48;

/// One cell with an edge on the intersection segment
#[derive(Clone, Copy, Debug)]
pub struct IntersectionCell {
    pub cell: usize,
    /// The fracture of the cell
    pub fracture: usize,
    /// Below 0: the flow goes out of the intersection into the cell
    pub outflow_product: f64,
    /// Square of the speed in the cell at the intersection point
    pub squared_speed: f64,
    /// The first fracture of the cell node which is not on the intersection
    pub third_node_fracture: usize,
}

impl IntersectionCell {
    fn is_outgoing(&self) -> bool {
        self.outflow_product < 0.0
    }

    fn speed(&self) -> f64 {
        self.squared_speed.sqrt()
    }
}

/// Complete mixing: an out-going cell, with a probability in proportion to
/// its speed (to the square of its speed for 4 or more out-going cells).
/// Without an out-going cell, the fastest cell other than `previous`.
/// Returns an index of `cells`.
pub fn complete_mixing(
    cells: &[IntersectionCell],
    previous: usize,
    random: &mut Drand48,
) -> Option<usize> {
    let outgoing: Vec<usize> = (0..cells.len())
        .filter(|&index| cells[index].is_outgoing())
        .collect();
    match outgoing.len() {
        0 => fastest_except(cells, previous),
        1 => Some(outgoing[0]),
        _ => Some(sample_outgoing(cells, &outgoing, random.next_f64())),
    }
}

/// The fastest cell other than `previous`, if a cell has a speed above 0
fn fastest_except(cells: &[IntersectionCell], previous: usize) -> Option<usize> {
    let mut fastest = None;
    let mut highest = 0.0;
    for (index, cell) in cells.iter().enumerate() {
        if cell.squared_speed > highest && index != previous {
            highest = cell.squared_speed;
            fastest = Some(index);
        }
    }
    fastest
}

/// The first out-going cell whose cumulative weight fraction is at least
/// the random number; the weight is the speed for up to 3 cells and the
/// square of the speed for more
fn sample_outgoing(cells: &[IntersectionCell], outgoing: &[usize], random_number: f64) -> usize {
    let weight = |index: usize| {
        if outgoing.len() <= 3 {
            cells[index].speed()
        } else {
            cells[index].squared_speed
        }
    };
    let total: f64 = outgoing.iter().map(|&index| weight(index)).sum();
    let mut cumulative = 0.0;
    for &index in &outgoing[..outgoing.len() - 1] {
        cumulative += weight(index);
        if random_number <= cumulative / total {
            return index;
        }
    }
    outgoing[outgoing.len() - 1]
}

/// Streamline routing. For two out-going cells: if one of them is in the
/// fracture of the incoming cell (continuous case), the particle goes into
/// the other one with the probability (out-going speed / incoming speed),
/// and always if that speed is the larger one; if both are in one fracture
/// (discontinuous case), the rule of Kang et al. In the other cases as
/// complete mixing. `previous_cell` is the cell of the particle before the
/// intersection, in `previous_fracture`.
pub fn streamline_routing(
    cells: &[IntersectionCell],
    previous: usize,
    previous_cell: usize,
    previous_fracture: usize,
    random: &mut Drand48,
) -> Option<usize> {
    let outgoing: Vec<usize> = (0..cells.len())
        .filter(|&index| cells[index].is_outgoing())
        .collect();
    let [first, second] = match outgoing[..] {
        [] => return fastest_except(cells, previous),
        [only] => return Some(only),
        [first, second] => [first, second],
        _ => return Some(sample_outgoing(cells, &outgoing, random.next_f64())),
    };

    // The opposite cell: an out-going cell in the fracture of the incoming
    // cell. The incoming speed: of the previous cell, else of a cell in
    // the fracture of the incoming cell which is not out-going.
    let mut opposite = None;
    let mut incoming_squared_speed = 0.0;
    for (index, cell) in cells.iter().enumerate() {
        let is_out = index == first || index == second;
        if cell.third_node_fracture == previous_fracture && is_out {
            opposite = Some(index);
        }
        if cell.third_node_fracture == previous_fracture && !is_out {
            incoming_squared_speed = cell.squared_speed;
        }
        if cell.cell == previous_cell {
            incoming_squared_speed = cell.squared_speed;
        }
    }
    let is_continuous = cells[first].fracture != cells[second].fracture;
    let random_number = random.next_f64();
    let incoming_speed = incoming_squared_speed.sqrt();

    if is_continuous {
        // DFNTrans: without an opposite cell, the second cell is the opposite one
        let (opposite, adjacent) = if opposite == Some(first) {
            (first, second)
        } else {
            (second, first)
        };
        if incoming_squared_speed <= cells[adjacent].squared_speed
            || random_number <= cells[adjacent].speed() / incoming_speed
        {
            Some(adjacent)
        } else {
            Some(opposite)
        }
    } else {
        let (speed_first, speed_second) = (cells[first].speed(), cells[second].speed());
        let (faster, slower) = if cells[second].squared_speed >= cells[first].squared_speed {
            (second, first)
        } else {
            (first, second)
        };
        if incoming_speed >= (speed_first + speed_second) / 2.0 {
            // On the faster incoming fracture
            if random_number >= cells[faster].speed() / incoming_speed {
                Some(slower)
            } else {
                Some(faster)
            }
        } else {
            // On the slower incoming fracture
            let opposite_speed = speed_first + speed_second - incoming_speed;
            let (larger, smaller) = if cells[second].squared_speed <= cells[first].squared_speed {
                (first, second)
            } else {
                (second, first)
            };
            if random_number <= (cells[larger].speed() - opposite_speed) / incoming_speed
                && cells[larger].speed() >= opposite_speed
            {
                Some(larger)
            } else {
                Some(smaller)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(index: usize, outflow_product: f64, squared_speed: f64) -> IntersectionCell {
        IntersectionCell {
            cell: index,
            fracture: index / 2,
            outflow_product,
            squared_speed,
            third_node_fracture: index / 2,
        }
    }

    #[test]
    fn complete_mixing_takes_the_only_outgoing_cell() {
        let cells = [cell(0, 1.0, 4.0), cell(1, -1.0, 1.0), cell(2, 1.0, 9.0)];
        let mut random = Drand48::new(1);
        assert_eq!(complete_mixing(&cells, 0, &mut random), Some(1));
    }

    #[test]
    fn complete_mixing_without_outgoing_cell_takes_the_fastest_other_cell() {
        let cells = [cell(0, 1.0, 9.0), cell(1, 1.0, 1.0), cell(2, 1.0, 4.0)];
        let mut random = Drand48::new(1);
        assert_eq!(complete_mixing(&cells, 0, &mut random), Some(2));
    }

    #[test]
    fn complete_mixing_samples_in_proportion_to_the_speed() {
        // Speeds 1 and 3: 3/4 of the particles into the second cell
        let cells = [cell(0, -1.0, 1.0), cell(1, -1.0, 9.0)];
        let mut random = Drand48::new(7);
        let count = 100_000;
        let second = (0..count)
            .filter(|_| complete_mixing(&cells, 0, &mut random) == Some(1))
            .count();
        let fraction = second as f64 / count as f64;
        assert!((fraction - 0.75).abs() < 0.01, "{fraction}");
    }

    #[test]
    fn streamline_routing_goes_into_a_faster_adjacent_branch() {
        // Incoming cell 0 (fracture 0), opposite cell 1 (fracture 0),
        // adjacent cell 2 (fracture 1) faster than the incoming flow
        let cells = [
            cell(0, 1.0, 1.0),
            cell(1, -1.0, 0.25),
            IntersectionCell {
                fracture: 1,
                third_node_fracture: 1,
                ..cell(2, -1.0, 4.0)
            },
        ];
        let mut random = Drand48::new(3);
        for _ in 0..100 {
            assert_eq!(streamline_routing(&cells, 0, 0, 0, &mut random), Some(2));
        }
    }
}
