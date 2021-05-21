import json
import eval7
import random

def calcualte_strength(hole, iters): 
        '''
        A Monte Carlo method meant to estimate the win probability of a pair of 
        hole cards. Simlulates 'iters' games and determines the win rates of our cards

        Arguments:
        hole: a list of our two hole cards
        iters: a integer that determines how many Monte Carlo samples to take
        '''

        deck = eval7.Deck() #eval7 object!
        hole_cards = [eval7.Card(card) for card in hole] #card objects, used to evaliate hands

        for card in hole_cards: #remove cards that we know about! they shouldn't come up in simulations
            deck.cards.remove(card)

        score = 0

        for _ in range(iters): #take 'iters' samples
            deck.shuffle() #make sure our samples are random

            _COMM = 5 #the number of cards we need to draw
            _OPP = 2

            draw = deck.peek(_COMM + _OPP)

            opp_hole = draw[: _OPP]
            community = draw[_OPP: ]

            our_hand = hole_cards + community #the two showdown hands
            opp_hand = opp_hole + community

            our_hand_value = eval7.evaluate(our_hand) #the ranks of our hands (only useful for comparisons)
            opp_hand_value = eval7.evaluate(opp_hand)

            if our_hand_value > opp_hand_value: #we win!
                score += 2
            
            elif our_hand_value == opp_hand_value: #we tie.
                score += 1
            
            else: #we lost....
                score += 0
        
        hand_strength = score / (2 * iters) #this is our win probability!

        return hand_strength

def allocation_generator():
    hand_percentages = {}
    ranks = ['2', '3', '4', '5', '6', '7', '8', '9', 'T', 'J', 'Q', 'K', 'A']
    for is_suit_same in ['s','d']: #s for same, d for different
        for rank1 in ranks:
            for rank2 in ranks:
                if is_suit_same == 's':
                    card1 = rank1 + 's'
                    card2 = rank2 + 's'
                else:
                    card1 = rank1 + 'h'
                    card2 = rank2 + 's'
                if card1 != card2:
                    hole_cand = (card1, card2)
                    hand_percentages[rank1+rank2+is_suit_same] = calcualte_strength(hole_cand, 50000)
    hand_percentages = {k: v for k, v in sorted(hand_percentages.items(), key=lambda item: item[1])}
    print(hand_percentages)
    with open('hand_percentages.txt', 'w') as outfile:
        json.dump(hand_percentages, outfile)

allocation_generator()