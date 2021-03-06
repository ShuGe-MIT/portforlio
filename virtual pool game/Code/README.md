**Project specification**

Introduction:
=======================================================================================
This project is a virtual pool game that allows two players to game together at different locations. Each of them simply needs an ESP32 to perform cue hits. They can also have a laptop to display the webpage, but this is optional, and it is just to see the animation of the balls moving after a shot.

Video Demo:
---------------------------------------------------------------------------------------

Video Demo:https://youtu.be/dQsGcX01YW8

*Background Music in Acknowledgments*

In this video, we first access the front end webpage to enter a gameid and usernames of the two players. This lets us view the pool table on the webpage.

Then, each of the players keys in their username and the gameid on the ESP32. They take turns hitting the ball -- when it is their turn, the pool table will load on the screen. Otherwise, there will be a message -- "Waiting for other player..." displayed. Notice that if a player fails to hit any ball in or pockets the cue ball, it will be the other player's turn. Otherwise, if they pocket a ball, the next turn is still theirs. Also, when the cue ball is pocketed, it reappears on the table. 

To hit the ball, the players rotate the ESP32 to change the angle of the cue stick and press the button to select the angle. After that, they hit the cue ball by mimicing a cue shot and pushing the ESP32 forward. The hit is then registered on the front end, and we see the animation displayed on the screen.

When the game ends, there is a pop-up on the webpage saying that the game ended, and both players' ESP32 restarts automatically.

Technical Overview:
=======================================================================================
System block diagram
---------------------------------------------------------------------------------------
<figure>
    <img src="https://i.ibb.co/Q6dwh8N/189139048-1775900755926274-4292546695778922304-n-1.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">
        block diagram
    </figcaption>
</figure>

During gameplay, after a hit, the ESP sends a POST request to the server, consisting of (username, gameid, stick_velocity, stick_angle). The server stores this tuplet in the database. 

The front end makes periodic GET requests to the server to check for new hits. Once a new hit is stored in the database, the front end can retrieve it and simulate the ball hit, and this will be animated on the webpage. Then, after the simulation is complete, the front end sends a POST request to the server consisting of a json string with information on the ball positions and who the next player is. The string is stored in the database. 

After the hit, the ESP makes periodic GET requests to the server to check if there is a new string in the database. It parses the json string to check if it is the next player -- if it is, it displays the balls on a mini pool table on the LCD screen, and allows the player to make the next hit.


Parts List
---------------------------------------------------------------------------------------
* ESP32
* IMU
* Buzzer
* LCD Screen
* Buttons


ESP32 Details:
=======================================================================================

Overall operations of ESP:
---------------------------------------------------------------------------------------

**Login stage**

1. Short press to start entering username.

2. Tilt the board to select alphabets/numbers.
   1. short press to confirm current character and move on with the next character.
   2. long press to enter the confirmation page.
3. confirmation page for username: long press to confirm and short press re-enter username.

4. Tilt the board to select alphabets/numbers.
   1. short press to confirm current character and move on with the next character.
   2. long press to enter the confirmation page.

5. confirmation page for game id: long press to confirm and short press re-enter game id.

**Game play stage**

6. if it is the player's turn:
   1. Pool table and balls are displayed on LCD. Cue ball is white, ball 8 is black, other balls are yellow.
   2. Rotate the board to adjust the angle of the cue stick.
   3. Short press to confirm the angle.
   4. Perform hit motion.
   5. If successful shot, continue with step 2 to perform next shot; otherwise, go to stage 7.

7. if it is not the player's turn:
   1. Nothing is displayed except for a message to wait for the next player.

**Restart**

8. Long press the right button to restart (go back to stage 1) at any point of the game.

9. Automatically returns to stage 1 when a game ends.

Operate ESP: https://youtu.be/IMU-J82otLo

Hardware Specification:
---------------------------------------------------------------------------------------

Hardware we have is an LCD screen and two buttons: the right one is only used for restarting during any time of the game when long-pressed; the left one is for all other button controls.

<figure>
    <img src="https://i.ibb.co/dLLL1Tm/IMG-2019.jpg"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">
        hard ware
    </figcaption>
</figure>



Comparison between pool game on LCD and the frontend website for the same state of the game.

<figure>
    <img src="https://i.ibb.co/kyzRB9S/IMG-2018.jpg"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">
        LCD pool table
    </figcaption>
</figure>

<figure>
    <img src="https://i.ibb.co/bXmkKN8/Screenshot-2021-05-20-172905.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">
        website pool table
    </figcaption>
</figure>


State machine for ESP
---------------------------------------------------------------------------------------
<figure>
    <img src="https://i.ibb.co/T4Wnhzn/Screenshot-2021-05-20-at-8-48-11-PM.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">
        State machine
    </figcaption>
</figure>


Arduino file explanation:
---------------------------------------------------------------------------------------

**cue_motion.ino**

This file contains the state machine for the ESP part of our project. The ESP starts off in the "Key username" state. Then, after keying in the username, a long button press changes the state to "Key gameid". In both of these states, users are able to input username and gameid through InputGetter.h. 

Another long button press confirms the gameid and moves on to the "Check current player and game state" state. In this state, periodic requests check_game_end() and retrieve() are sent to the server to check if the game has ended or if it is the player's turn to hit the ball.

If the game has ended, it sends the player back to the "Key username" state, essentially "restarting" the ESP. If it is the player's turn to hit the ball, it moves on to the "Adjust angle of cue stick" state. In this state, the IMU detects the angular acceleration of the board using Gyroscope. If angular acceleration is greater than a threshold, we increase the global variable ???angle??? by step_ang = 0.1 radian; if the acceleration is smaller than a threshold, we decrease the ???angle??? by 0.1. Therefore, the cue will rotate in the same direction of rotation of the ESP.

Once the user fixes the angle, the state machine goes into the hit cue ball stage. The accelerometer in IMU detects acceleration. Hit motion is detected if acceleration is greater than 1g. The buzzer will play some notes by calling run_buzzer function as discussed below. If acceleration is less than or equal to 1g, no hit motion is detected. The user can give another jerk to the board until a hit is detected or short press to re-adjust angles.

If a hit is detected, we will post the cue velocity (proportional to acceleration with scale of 1), cue angle from horizontal, game id, and user name to the server. The state machine goes back to ???check current player and game state??? stage.


**InputGetter.h**

InputGetter class is used to collect user inputs for username and game id. It is similar to the WikipediaGetter in exercise 5. The accelerometer detects tilting of the board and allows switching between characters. It keeps track of the button inputs through a state variable input_state. It updates global variables username and gameid once the user confirms the inputs, which are later sent to the server in POST requests.

**Button.h**

Button class is used to detect long presses and short presses.

**buzzer.ino**

The run_buzzer() function takes the magnitude of hit force applied as an input and creates a buzzer sound with the number of notes played proportional to the magnitude of the force.

**server_communication.ino**

There are three functions in this file:

post_hit() sends POST request to server file hitball.py after a hit motion is detected. The request contains user name (username), game id (game_id), initial velocity of the cue ball proportional to acceleration of the hit motion (stick_vel), and angle of the stick from horizontal (stick_ang).

check_game_end() sends GET request to the server file result_post.py periodically to continuously check if the game has ended. It gets a Boolean value from the server, True for game end.

retrieve() sends a GET request to the server to get current player, positions of cue ball, ball 8 (black ball), and other balls that are still on the pool table. If the user is the current player, retrieve function will draw pool table and all available balls on the LCD. However, if the user is not the current player, the LCD will only display waiting message, without the pool table.

Server Details:
=======================================================================================

The server handles both POST and GET requests from both the ESP and the front end. It serves as an intermediary between the front end and the ESP as well as maintains a database so that both old games and score records can be stored. 

Providing Front End with Memory
---------------------------------------------------------------------------------------

After the users log in from the website, the front end sends a POST request to server containing gameid and usernames of the two players. The server checks player_table (see below) to see if there's any entry with the same gameid. If so, it means it's an old game, and server returns latest stored ball positions from game_table (explained later); if not, it means it's a new game, and server inserts a new entry containing gameid and usernames into player_table.

 Game ID | Player 1 | Player 2
:-------:|:--------:|:--------:
  111111 | tim      | mit
  [player_table]

Storing Hit Information from ESP
---------------------------------------------------------------------------------------

After every hit motion on the ESP, the ESP sends a POST request to the server containing gameid, player name, stick angle and velocity; the server file checks whether it is this player's turn by getting the latest entry from next_player_table (explained later). If it is this player's turn, server file inserts a new entry into hit_table (see below).

 Game ID | Player Name | Stick Angle | Stick Velocity
:-------:|:-----------:|:-----------:|:--------:
  111111 |   timbeaver |  3          | 6
  [hit_table]

Providing Front End with Hit Information
---------------------------------------------------------------------------------------

The front end sends periodic GET requests to server to check if there is a new entry being inserted into hit_table. If a new entry is inserted, front end will use the stick angle and stick velocity to run the simulation.

Storing Updated Game State from Front End
---------------------------------------------------------------------------------------

After player makes a hit and front end finishes running the simulation, front ends sends a POST request containing information on all the ball positions and the identity of the next player to server. Server stores the updated ball positions in game_table (see below) and the username of the next player in next_player_table (see below).

 Game ID | Game State                                      | Timing
:-------:|:-----------------------------------------------:|:----------------------:
  111111 | json object with ball positions and next player | 2021-5-20 08:23:19.120
  [game_table]

 Game ID | Next Player | Timing
:-------:|:-----------:|:----------------------:
  111111 | timbeaver   | 2021-5-20 08:23:19.120
  [next_player_table]

Providing ESP with Ball Positions
---------------------------------------------------------------------------------------

After each hit, ESP sends a GET request to server containing gameid. Server retrieves latest entry from game_table under the same gameid, parses the json object and returns the ball positions to ESP.

Storing Game Result from Front End
---------------------------------------------------------------------------------------

After the game ends, the front end sends a POST request to server containing gameid, winner and loser. Server stores these information in results_table (see below).

 Game ID | Winner      | Loser
:-------:|:-----------:|:----------------------:
  111111 | mit         | tim
  [results_table]

Telling ESP Whether Game Has Ended
---------------------------------------------------------------------------------------

Right after log in and after every hit, ESP sends a GET request to server containing gameid. Server checks whether there's an entry in results_table under the same gameid. If so, it means game has already ended, and server returns True; if not, it means game is ongoing, and server returns False.

Providing Front End with Score Board Information
---------------------------------------------------------------------------------------

After players log in from website, front end sends a POST request containing usernames of the two players to server. Server checks the results_table to see how many games each player has played previously and how many games each has won. Server then returns this information to front end for it to calculate and display win rates of the players.

Front end:
=======================================================================================
The below link will take you to the frontend website of our pool game project.

[Pool Game Frontend Website](https://pool-game-project.herokuapp.com/)

*To enable the game, add the website to the allowed website list of insecure content settings in Chrome. This will allow the page to make XMLHttpRequests to the server.*

Our frontend is built with Javascript, HTML, and CSS. We deplyed it using heroku so that it can run on any browsers without any installations.

We devide our frontend website to the following parts:
- Log in
- Pool Table
- Hit Message
- Game End Pop-up

Log in
---------------------------------------------------------------------------------------
The users are asked to log in to a game once they open up the website. 
By inputting a 6-digit game ID and the name of two players, they can either initiate a new game or recover an unfinished old game.
This is achieved by sending a POST request to the server to search in the game database.
If there exist an old game in the database, game state information regarding the current player and
the position of balls will be returned and displayed on the pool table.

<figure>
    <img src="https://i.ibb.co/1Rpb5sD/Screen-Shot-2021-05-20-at-11-22-57-AM.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">Log in Page</figcaption>
</figure>

Pool Table
---------------------------------------------------------------------------------------
The pool table is where the movement of balls are animated and displayed.
When the current player makes a hit on the ESP 32, 
the angle and acceleration of the hit motion is sent to the server, 
to which frontend makes a periodic GET request to retrieve the latest hit.
After each hit, the updated game state information of ball positions and current player are sent to be stored at the database.

The players' scoreboard is shown at the bottom left corner, displaying the game history of players. 
At the beggining of the game, this information is retrived from the players' record database on the server.
At the end of the game, the winner and loser of the game will be POSTed to the database to store.
For new players, the win rate is set to be "???".
<figure>
    <img src="https://i.ibb.co/VJ3dQpq/Screen-Shot-2021-05-20-at-11-23-15-AM.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">Pool Table Page</figcaption>
</figure>

Hit Message
---------------------------------------------------------------------------------------
Hit messages are displayed when there's a target ball being hit into the pocket 
and when a player hit balls into pockets in a roll.
<figure>
    <img src="https://i.ibb.co/VC9Kp8t/Screen-Shot-2021-05-20-at-11-23-32-AM.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">Hit Message Display</figcaption>
</figure>

Game End Pop-up
---------------------------------------------------------------------------------------
When the game ends, a pop-up will annouce the winner of the game and ask the player if they want to restart a new game.
If players confirm, the website will redirect to the game log in page.
If players deny, the website will close automatically.
For both options, a second pop-up will appear to ask players to confirm their selection.
<figure>
    <img src="https://i.ibb.co/5kfCW5S/Screen-Shot-2021-05-20-at-11-24-16-AM.png"
     style="float: left; margin-right: 10px;" />
    <figcaption style="text-align: center;font-weight:bold;">Game End Pop-up</figcaption>
</figure>

Acknowledgments
=======================================================================================
We want to thank our 6.08 instructors, Karl, Varnika, and Annika for their guidance and support.

Background Music for Video Demo:
Bliss by Luke Bergs | https://soundcloud.com/bergscloud/
Creative Commons - Attribution-ShareAlike 3.0 Unported
https://creativecommons.org/licenses/by-sa/3.0/
Music promoted by https://www.chosic.com/
