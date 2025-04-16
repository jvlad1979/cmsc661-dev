def calculate_price_per_oz_alcohol(item_volume, abv_percent, cost, quantity, unit_type):
    """
    Calculates the price per ounce of pure alcohol in a beverage,
    handling different unit types (ml, oz) and quantity.

    Args:
        item_volume (float): The volume of a single item (bottle, can) in ml or oz.
        abv_percent (float): The alcohol by volume (ABV) as a percentage (e.g., 40 for 40%).
        cost (float): The cost of the beverage.
        quantity (int): The number of items (bottles, cans).  For wine/liquor, this will typically be 1.
        unit_type (str): A string indicating the unit type:
            'ml' for milliliters,
            'oz' for fluid ounces.

    Returns:
        float: The price per ounce of pure alcohol, or None if the input is invalid.
    """
    if not all(isinstance(arg, (int, float)) and arg >= 0 for arg in (item_volume, abv_percent, cost)) or not isinstance(quantity, int) or quantity <= 0:
        return None  # Return None for invalid input

    ml_to_oz = 0.033814  # Conversion factor

    if unit_type not in ('ml', 'oz'):
        return None  # Return None for invalid unit type

    # Calculate the volume of alcohol per item
    if unit_type == 'ml':
        alcohol_volume_ml_per_item = item_volume * (abv_percent / 100)
        alcohol_volume_oz_per_item = alcohol_volume_ml_per_item * ml_to_oz
    elif unit_type == 'oz':
        alcohol_volume_oz_per_item = item_volume * (abv_percent / 100)
    else:
        return None

    # Calculate the total volume of alcohol
    total_alcohol_volume_oz = alcohol_volume_oz_per_item * quantity

    if total_alcohol_volume_oz == 0:
        return None  # Return None if the alcohol volume is 0

    price_per_oz_alcohol = cost / total_alcohol_volume_oz
    return price_per_oz_alcohol


def main():
    """
    Main function to run the calculator with user inputs.
    """
    print("Welcome to the Alcohol Price Calculator!")
    print("Enter the beverage details to calculate the price per ounce of pure alcohol.")

    while True:
        unit_type = input("Enter the unit type ('ml', 'oz'): ").lower()
        if unit_type not in ('ml', 'oz'):
            print("Invalid unit type. Please enter 'ml' or 'oz'.")
            continue

        if unit_type == 'ml':
            beverage_type = "wine/liquor"
        else:
            beverage_type = "beer/cider/etc."

        break

    while True:
        try:
            item_volume = float(input(f"Enter the volume of a single {beverage_type} item in {unit_type}: "))
            if item_volume < 0:
                print("Volume must be a non-negative number. Please try again.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter a valid number for volume.")

    while True:
        try:
            if beverage_type == "wine/liquor":
                quantity = 1
                print("Quantity is set to 1 as it is a single bottle.")
            else:
                quantity = int(input("Enter the quantity of items (cans/bottles) in the pack: "))
                if quantity < 1:
                    print("Quantity must be a positive integer. Please try again.")
                    continue
            break
        except ValueError:
            print("Invalid input. Please enter a valid number for quantity.")

    while True:
        try:
            abv_percent = float(input("Enter the ABV percentage (e.g., 40 for 40%): "))
            if abv_percent < 0:
                print("ABV must be a non-negative number. Please try again.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter a valid number for ABV.")

    while True:
        try:
            cost = float(input(f"Enter the cost of the {beverage_type}: "))
            if cost < 0:
                print("Cost must be a non-negative number. Please try again.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter a valid number for cost.")

    result = calculate_price_per_oz_alcohol(item_volume, abv_percent, cost, quantity, unit_type)

    if result is None:
        print("Invalid input: Please ensure all values are non-negative, and that the volume and ABV are not zero, and the unit type is valid.")
    else:
        print(f"The price per ounce of pure alcohol is: ${result:.2f}")

    print("Thank you for using the Alcohol Price Calculator!")


if __name__ == "__main__":
    main()
