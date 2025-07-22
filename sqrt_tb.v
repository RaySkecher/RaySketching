`timescale 1ns/1ps

/*
 * Module: sqrt_tb
 *
 * Description:
 *   A testbench for the 16-bit integer square root calculator (sqrt).
 */
module sqrt_tb;

// Testbench signals
reg clk;
reg rst_n;
reg start;
reg [15:0] data_in;

// DUT outputs
wire [7:0] data_out;
wire done;

// Instantiate the Unit Under Test (UUT)
sqrt uut (
    .clk(clk),
    .rst_n(rst_n),
    .start(start),
    .data_in(data_in),
    .data_out(data_out),
    .done(done)
);

// Clock generation (100MHz)
initial begin
    clk = 0;
    forever #5 clk = ~clk;
end

// Test sequence
initial begin
    $display("Starting sqrt testbench...");

    // Initialize signals and apply reset
    start <= 0;
    data_in <= 0;
    rst_n <= 0;
    #20;
    rst_n <= 1;
    #10;

    // --- Test Cases ---

    // Test case 1: sqrt(0) -> 0
    test_value(16'd0, 8'd0);

    // Test case 2: sqrt(64) -> 8
    test_value(16'd64, 8'd8);

    // Test case 3: sqrt(4096) -> 64
    test_value(16'd4096, 8'd64);

    // Test case 4: sqrt(5000) -> 70 (floor of 70.71)
    test_value(16'd5000, 8'd70);

    // Test case 5: sqrt(65535) -> 255 (floor of 255.99)
    test_value(16'd65535, 8'd255);

    // Test case 6: sqrt(2) -> 1
    test_value(16'd2, 8'd1);

    // Test case 7: A value in between -> sqrt(1000) -> 31
    test_value(16'd1000, 8'd31);

    $display("All test cases passed!");
    $finish;
end

// Task to apply a value and check the result
task test_value;
    input [15:0] input_val;
    input [7:0]  expected_val;
    begin
        @(posedge clk);
        data_in <= input_val;
        start <= 1;
        $display("Testing sqrt(%0d)...", input_val);
        @(posedge clk);
        start <= 0;

        // Wait for the done signal to be asserted
        while (!done) begin
            @(posedge clk);
        end

        // Check the result on the same cycle 'done' is high
        if (data_out !== expected_val) begin
            $display("TEST FAILED: Input=%0d, Expected=%0d, Got=%0d", input_val, expected_val, data_out);
            $finish;
        end else begin
            $display("Test Passed: Input=%0d, Sqrt=%0d", input_val, data_out);
        end
        
        // Wait for done to go low before starting next test
        @(posedge clk);
    end
endtask

endmodule 