/*
 * Module: sqrt
 *
 * Description:
 *   A 16-bit integer square root calculator. It uses a sequential, restoring
 *   shift-and-subtract algorithm.
 *
 *   The calculation takes 8 clock cycles to complete.
 *
 * Interface:
 *   clk     - Input clock
 *   rst_n   - Asynchronous active-low reset
 *   start   - Input signal to begin the square root calculation.
 *             The module must be idle (done=1 or initial state) to start.
 *   data_in - 16-bit input value (radicand).
 *   data_out- 8-bit output for the integer part of the square root.
 *   done    - Output signal that goes high for one clock cycle when the
 *             calculation is complete.
 */
module sqrt(
    input wire clk,
    input wire rst_n,
    input wire start,
    input wire [15:0] data_in,
    output reg [7:0] data_out,
    output reg done
);

// State machine states
localparam STATE_IDLE = 2'b00;
localparam STATE_CALC = 2'b01;
localparam STATE_DONE = 2'b10;

reg [1:0] state;

// Internal registers for calculation
reg [17:0] remainder;
reg [15:0] radicand_shifted;
reg [7:0]  current_root;
reg [3:0]  iter_count; // 8 iterations for 8-bit result

// Temporary variables for calculation logic
reg [17:0] dividend;
reg [9:0]  trial_subtrahend;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        // Asynchronous reset
        state <= STATE_IDLE;
        data_out <= 8'h00;
        done <= 1'b0;
        remainder <= 18'h0;
        radicand_shifted <= 16'h0000;
        current_root <= 8'h00;
        iter_count <= 4'd0;
    end else begin
        // Default assignments to avoid latches in combinational paths
        done <= 1'b0;

        case (state)
            STATE_IDLE: begin
                if (start) begin
                    // Latch input and initialize for calculation
                    radicand_shifted <= data_in;
                    remainder <= 18'h0;
                    current_root <= 8'h00;
                    iter_count <= 4'd8; // Countdown from 8
                    state <= STATE_CALC;
                end
            end

            STATE_CALC: begin
                if (iter_count > 0) begin
                    // Perform one iteration of the restoring square root algorithm
                    
                    // The next part of the number to consider
                    dividend = {remainder[15:0], radicand_shifted[15:14]};
                    
                    // The value to attempt to subtract: (root << 2) | 1
                    trial_subtrahend = {current_root, 2'b01};

                    if (dividend >= trial_subtrahend) begin
                        remainder <= dividend - trial_subtrahend;
                        current_root <= {current_root[6:0], 1'b1};
                    end else begin
                        remainder <= dividend; // No change in remainder if subtraction doesn't happen
                        current_root <= {current_root[6:0], 1'b0};
                    end

                    // Shift radicand for next iteration
                    radicand_shifted <= radicand_shifted << 2;
                    iter_count <= iter_count - 1;
                end

                if (iter_count == 0) begin // Transition to DONE after last iteration
                    state <= STATE_DONE;
                end
            end

            STATE_DONE: begin
                data_out <= current_root;
                done <= 1'b1;
                state <= STATE_IDLE;
            end

            default: begin
                state <= STATE_IDLE;
            end
        endcase
    end
end

endmodule 