<?php
// MJPEG stream — pushes capture.jpg to the browser as a continuous
// multipart stream. The browser opens one connection and frames
// arrive as fast as they are written. No polling, no repeated requests.
//
// Apache timeout must be extended (see efinder.conf) otherwise Apache
// will kill the connection after the default TimeOut (usually 60s).

$image = '/dev/shm/efinder_live.jpg';
$interval_us = 500000; // 0.5 seconds between frames — tune as needed

header('Content-Type: multipart/x-mixed-replace; boundary=frame');
header('Cache-Control: no-store');
header('Connection: close');

// Disable PHP output buffering so frames reach the browser immediately.
while (ob_get_level() > 0) {
    ob_end_clean();
}

while (true) {
    if (file_exists($image)) {
        $data = @file_get_contents($image);
        if ($data !== false && strlen($data) > 0) {
            echo "--frame\r\n";
            echo "Content-Type: image/jpeg\r\n";
            echo "Content-Length: " . strlen($data) . "\r\n\r\n";
            echo $data . "\r\n";
            flush();
        }
    }
    usleep($interval_us);

    // Exit cleanly if the client has disconnected.
    if (connection_aborted()) {
        break;
    }
}
