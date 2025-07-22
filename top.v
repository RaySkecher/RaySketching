/*
 * Module: top
 *
 * Description:
 *   A top-level wrapper for the 16-bit integer square root calculator (sqrt).
 *   This module connects the sqrt core to the top-level I/O of an FPGA,
 *   specifically targeting the LEDs on a board like the Basys 3. This ensures
 *   that the synthesis tool (e.g., Vivado) does not trim the logic away.
 *
 *   This is intended to be the top module for synthesis.
 *
 *   - sw[15:0] can be connected to the 16 switches.
 *   - led[7:0] will display the 8-bit result of the square root.
 *   - led[8] will light up when the calculation is done.
 *
 */
module top(
    input         clk,
    input         rst,
    input         start,
    input  [15:0] sw,
    output [15:0] led
);

    wire [7:0]  sqrt_result;
    wire        sqrt_done;
    wire        rst_n;

    // Invert the active-high reset from the button for the active-low reset of the core
    assign rst_n = ~rst;

    // Instantiate the sqrt module
    sqrt sqrt_inst (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .data_in(sw),
        .data_out(sqrt_result),
        .done(sqrt_done)
    );

    // Assign outputs to LEDs
    // Lower 8 LEDs show the square root result
    assign led[7:0]   = sqrt_result;
    // LED[8] shows the done signal
    assign led[8]     = sqrt_done;
    // The rest of the LEDs are off
    assign led[15:9]  = 7'h0;

endmodule 